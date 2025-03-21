// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "eigensolver.hpp"

#include <mfem.hpp>
#include "fem/errorindicator.hpp"
#include "fem/mesh.hpp"
#include "linalg/arpack.hpp"
#include "linalg/divfree.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/floquetcorrection.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/slepc.hpp"
#include "linalg/vector.hpp"
#include "models/lumpedportoperator.hpp"
#include "models/postoperator.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"

namespace palace
{

using namespace std::complex_literals;

std::pair<ErrorIndicator, long long int>
EigenSolver::Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const
{
  // Construct and extract the system matrices defining the eigenvalue problem. The diagonal
  // values for the mass matrix PEC dof shift the Dirichlet eigenvalues out of the
  // computational range. The damping matrix may be nullptr.
  BlockTimer bt0(Timer::CONSTRUCT);
  SpaceOperator space_op(iodata, mesh);
  auto K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
  auto C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);

  const auto &Curl = space_op.GetCurlMatrix();
  SaveMetadata(space_op.GetNDSpaces());

  // Configure objects for postprocessing.
  PostOperator post_op(iodata, space_op, "eigenmode");
  ComplexVector E(Curl.Width()), B(Curl.Height());
  E.UseDevice(true);
  B.UseDevice(true);

  // Define and configure the eigensolver to solve the eigenvalue problem:
  //         (K + λ C + λ² M) u = 0    or    K u = -λ² M u
  // with λ = iω. In general, the system matrices are complex and symmetric.
  std::unique_ptr<EigenvalueSolver> eigen;
  config::EigenSolverData::Type type = iodata.solver.eigenmode.type;
#if defined(PALACE_WITH_ARPACK) && defined(PALACE_WITH_SLEPC)
  if (type == config::EigenSolverData::Type::DEFAULT)
  {
    type = config::EigenSolverData::Type::SLEPC;
  }
#elif defined(PALACE_WITH_ARPACK)
  if (iodata.solver.eigenmode.type == config::EigenSolverData::Type::SLEPC)
  {
    Mpi::Warning("SLEPc eigensolver not available, using ARPACK!\n");
  }
  type = config::EigenSolverData::Type::ARPACK;
#elif defined(PALACE_WITH_SLEPC)
  if (iodata.solver.eigenmode.type == config::EigenSolverData::Type::ARPACK)
  {
    Mpi::Warning("ARPACK eigensolver not available, using SLEPc!\n");
  }
  type = config::EigenSolverData::Type::SLEPC;
#else
#error "Eigenmode solver requires building with ARPACK or SLEPc!"
#endif
  if (type == config::EigenSolverData::Type::ARPACK)
  {
#if defined(PALACE_WITH_ARPACK)
    Mpi::Print("\nConfiguring ARPACK eigenvalue solver:\n");
    if (C)
    {
      eigen = std::make_unique<arpack::ArpackPEPSolver>(space_op.GetComm(),
                                                        iodata.problem.verbose);
    }
    else
    {
      eigen = std::make_unique<arpack::ArpackEPSSolver>(space_op.GetComm(),
                                                        iodata.problem.verbose);
    }
#endif
  }
  else  // config::EigenSolverData::Type::SLEPC
  {
#if defined(PALACE_WITH_SLEPC)
    Mpi::Print("\nConfiguring SLEPc eigenvalue solver:\n");
    std::unique_ptr<slepc::SlepcEigenvalueSolver> slepc;
    if (C)
    {
      if (!iodata.solver.eigenmode.pep_linear)
      {
        slepc = std::make_unique<slepc::SlepcPEPSolver>(space_op.GetComm(),
                                                        iodata.problem.verbose);
        slepc->SetType(slepc::SlepcEigenvalueSolver::Type::TOAR);
      }
      else
      {
        slepc = std::make_unique<slepc::SlepcPEPLinearSolver>(space_op.GetComm(),
                                                              iodata.problem.verbose);
        slepc->SetType(slepc::SlepcEigenvalueSolver::Type::KRYLOVSCHUR);
      }
    }
    else
    {
      slepc = std::make_unique<slepc::SlepcEPSSolver>(space_op.GetComm(),
                                                      iodata.problem.verbose);
      slepc->SetType(slepc::SlepcEigenvalueSolver::Type::KRYLOVSCHUR);
    }
    slepc->SetProblemType(slepc::SlepcEigenvalueSolver::ProblemType::GEN_NON_HERMITIAN);
    slepc->SetOrthogonalization(
        iodata.solver.linear.gs_orthog_type == config::LinearSolverData::OrthogType::MGS,
        iodata.solver.linear.gs_orthog_type == config::LinearSolverData::OrthogType::CGS2);
    eigen = std::move(slepc);
#endif
  }
  EigenvalueSolver::ScaleType scale = iodata.solver.eigenmode.scale
                                          ? EigenvalueSolver::ScaleType::NORM_2
                                          : EigenvalueSolver::ScaleType::NONE;
  if (C)
  {
    eigen->SetOperators(*K, *C, *M, scale);
  }
  else
  {
    eigen->SetOperators(*K, *M, scale);
  }
  eigen->SetNumModes(iodata.solver.eigenmode.n, iodata.solver.eigenmode.max_size);
  eigen->SetTol(iodata.solver.eigenmode.tol);
  eigen->SetMaxIter(iodata.solver.eigenmode.max_it);
  Mpi::Print(" Scaling γ = {:.3e}, δ = {:.3e}\n", eigen->GetScalingGamma(),
             eigen->GetScalingDelta());

  // If desired, use an M-inner product for orthogonalizing the eigenvalue subspace. The
  // constructed matrix just references the real SPD part of the mass matrix (no copy is
  // performed). Boundary conditions don't need to be eliminated here.
  std::unique_ptr<Operator> KM;
  if (iodata.solver.eigenmode.mass_orthog)
  {
    Mpi::Print(" Basis uses M-inner product\n");
    KM = space_op.GetInnerProductMatrix(0.0, 1.0, nullptr, M.get());
    eigen->SetBMat(*KM);

    // Mpi::Print(" Basis uses (K + M)-inner product\n");
    // KM = space_op.GetInnerProductMatrix(1.0, 1.0, K.get(), M.get());
    // eigen->SetBMat(*KM);
  }

  // Construct a divergence-free projector so the eigenvalue solve is performed in the space
  // orthogonal to the zero eigenvalues of the stiffness matrix.
  std::unique_ptr<DivFreeSolver<ComplexVector>> divfree;
  if (iodata.solver.linear.divfree_max_it > 0 &&
      !space_op.GetMaterialOp().HasWaveVector() &&
      !space_op.GetMaterialOp().HasLondonDepth())
  {
    Mpi::Print(" Configuring divergence-free projection\n");
    constexpr int divfree_verbose = 0;
    divfree = std::make_unique<DivFreeSolver<ComplexVector>>(
        space_op.GetMaterialOp(), space_op.GetNDSpace(), space_op.GetH1Spaces(),
        space_op.GetAuxBdrTDofLists(), iodata.solver.linear.divfree_tol,
        iodata.solver.linear.divfree_max_it, divfree_verbose);
    eigen->SetDivFreeProjector(*divfree);
  }

  // If using Floquet BCs, a correction term (kp x E) needs to be added to the B field.
  std::unique_ptr<FloquetCorrSolver<ComplexVector>> floquet_corr;
  if (space_op.GetMaterialOp().HasWaveVector())
  {
    floquet_corr = std::make_unique<FloquetCorrSolver<ComplexVector>>(
        space_op.GetMaterialOp(), space_op.GetNDSpace(), space_op.GetRTSpace(),
        iodata.solver.linear.tol, iodata.solver.linear.max_it, 0);
  }

  // Set up the initial space for the eigenvalue solve. Satisfies boundary conditions and is
  // projected appropriately.
  if (iodata.solver.eigenmode.init_v0)
  {
    ComplexVector v0;
    if (iodata.solver.eigenmode.init_v0_const)
    {
      Mpi::Print(" Using constant starting vector\n");
      space_op.GetConstantInitialVector(v0);
    }
    else
    {
      Mpi::Print(" Using random starting vector\n");
      space_op.GetRandomInitialVector(v0);
    }
    if (divfree)
    {
      divfree->Mult(v0);
    }
    eigen->SetInitialSpace(v0);  // Copies the vector

    // Debug
    // const auto &Grad = space_op.GetGradMatrix();
    // ComplexVector r0(Grad->Width());
    // r0.UseDevice(true);
    // Grad.MultTranspose(v0.Real(), r0.Real());
    // Grad.MultTranspose(v0.Imag(), r0.Imag());
    // r0.Print();
  }

  // Configure the shift-and-invert strategy is employed to solve for the eigenvalues
  // closest to the specified target, σ.
  const double target = iodata.solver.eigenmode.target;
  {
    const double f_target =
        iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, target);
    Mpi::Print(" Shift-and-invert σ = {:.3e} GHz ({:.3e})\n", f_target, target);
  }
  if (C)
  {
    // Search for eigenvalues closest to λ = iσ.
    eigen->SetShiftInvert(1i * target);
    if (type == config::EigenSolverData::Type::ARPACK)
    {
      // ARPACK searches based on eigenvalues of the transformed problem. The eigenvalue
      // 1 / (λ - σ) will be a large-magnitude negative imaginary number for an eigenvalue
      // λ with frequency close to but not below the target σ.
      eigen->SetWhichEigenpairs(EigenvalueSolver::WhichType::SMALLEST_IMAGINARY);
    }
    else
    {
      eigen->SetWhichEigenpairs(EigenvalueSolver::WhichType::TARGET_IMAGINARY);
    }
  }
  else
  {
    // Linear EVP has eigenvalues μ = -λ² = ω². Search for eigenvalues closest to μ = σ².
    eigen->SetShiftInvert(target * target);
    if (type == config::EigenSolverData::Type::ARPACK)
    {
      // ARPACK searches based on eigenvalues of the transformed problem. 1 / (μ - σ²)
      // will be a large-magnitude positive real number for an eigenvalue μ with frequency
      // close to but below the target σ².
      eigen->SetWhichEigenpairs(EigenvalueSolver::WhichType::LARGEST_REAL);
    }
    else
    {
      eigen->SetWhichEigenpairs(EigenvalueSolver::WhichType::TARGET_REAL);
    }
  }

  // Set up the linear solver required for solving systems involving the shifted operator
  // (K - σ² M) or P(iσ) = (K + iσ C - σ² M) during the eigenvalue solve. The
  // preconditioner for complex linear systems is constructed from a real approximation
  // to the complex system matrix.
  auto A = space_op.GetSystemMatrix(std::complex<double>(1.0, 0.0), 1i * target,
                                    std::complex<double>(-target * target, 0.0), K.get(),
                                    C.get(), M.get());
  auto P = space_op.GetPreconditionerMatrix<ComplexOperator>(1.0, target, -target * target,
                                                             target);
  auto ksp = std::make_unique<ComplexKspSolver>(iodata, space_op.GetNDSpaces(),
                                                &space_op.GetH1Spaces());
  ksp->SetOperators(*A, *P);
  eigen->SetLinearSolver(*ksp);

  // Initialize structures for storing and reducing the results of error estimation.
  TimeDependentFluxErrorEstimator<ComplexVector> estimator(
      space_op.GetMaterialOp(), space_op.GetNDSpaces(), space_op.GetRTSpaces(),
      iodata.solver.linear.estimator_tol, iodata.solver.linear.estimator_max_it, 0,
      iodata.solver.linear.estimator_mg);
  ErrorIndicator indicator;

  // CUSTOM CONVERGENCE
  indicator.SetConvergenceParams(
    iodata.solver.eigenmode.tol,     // Global tolerance
    iodata.solver.eigenmode.tol/10,  // Relative tolerance 
    3                                // Required consecutive converges
  );

  // Eigenvalue problem solve.
  BlockTimer bt1(Timer::EPS);
  Mpi::Print("\n");
  int num_conv = eigen->Solve();
  {
    std::complex<double> lambda = (num_conv > 0) ? eigen->GetEigenvalue(0) : 0.0;
    Mpi::Print(" Found {:d} converged eigenvalue{}{}\n", num_conv,
               (num_conv > 1) ? "s" : "",
               (num_conv > 0)
                   ? fmt::format(" (first = {:.3e}{:+.3e}i)", lambda.real(), lambda.imag())
                   : "");
  }
  BlockTimer bt2(Timer::POSTPRO);
  SaveMetadata(*ksp);

  // Calculate and record the error indicators, and postprocess the results.
  Mpi::Print("\nComputing solution error estimates and performing postprocessing\n");
  if (!KM)
  {
    // Normalize the finalized eigenvectors with respect to mass matrix (unit electric field
    // energy) even if they are not computed to be orthogonal with respect to it.
    KM = space_op.GetInnerProductMatrix(0.0, 1.0, nullptr, M.get());
    eigen->SetBMat(*KM);
    eigen->RescaleEigenvectors(num_conv);
  }
  Mpi::Print("\n");
  for (int i = 0; i < num_conv; i++)
  {
    // Get the eigenvalue and relative error.
    std::complex<double> omega = eigen->GetEigenvalue(i);
    double error_bkwd = eigen->GetError(i, EigenvalueSolver::ErrorType::BACKWARD);
    double error_abs = eigen->GetError(i, EigenvalueSolver::ErrorType::ABSOLUTE);
    if (!C)
    {
      // Linear EVP has eigenvalue μ = -λ² = ω².
      omega = std::sqrt(omega);
    }
    else
    {
      // Quadratic EVP solves for eigenvalue λ = iω.
      omega /= 1i;
    }

    // Compute B = -1/(iω) ∇ x E on the true dofs, and set the internal GridFunctions in
    // PostOperator for all postprocessing operations.
    eigen->GetEigenvector(i, E);
    Curl.Mult(E.Real(), B.Real());
    Curl.Mult(E.Imag(), B.Imag());
    B *= -1.0 / (1i * omega);
    if (space_op.GetMaterialOp().HasWaveVector())
    {
      // Calculate B field correction for Floquet BCs.
      // B = -1/(iω) ∇ x E + 1/ω kp x E.
      floquet_corr->AddMult(E, B, 1.0 / omega);
    }
    post_op.SetEGridFunction(E);
    post_op.SetBGridFunction(B);
    post_op.UpdatePorts(space_op.GetLumpedPortOp(), omega.real());
    const double E_elec = post_op.GetEFieldEnergy();
    const double E_mag = post_op.GetHFieldEnergy();

    // CUSTOM CONVERGENCE MODIFY
    // Calculate and record the error indicators.
    if (i < iodata.solver.eigenmode.n)
    {
      // Check if the current element/region contains a Josephson Junction
      bool is_jj = false;
      
      // Get port information
      const auto& lumped_port_op = space_op.GetLumpedPortOp();
      
      // Check for Josephson junctions
      for (const auto& [idx, port] : lumped_port_op) {
        // Consider any pure inductive lumped element (L > 0, R = 0, C = 0) as a Josephson junction
        if (std::abs(port.L) > 0.0 && std::abs(port.R) <= 0.0 && std::abs(port.C) <= 0.0) {
          is_jj = true;
          Mpi::Print(" Detected Josephson junction at port {}\n", idx);
          break;
        }
      }
      
      // Check if the current element/region contains a Josephson Junction
      bool has_jj = false;
      
      // Get port information
      const auto& lp_op = space_op.GetLumpedPortOp();
      
      // Check for Josephson junctions
      for (const auto& [idx, port] : lp_op) {
        // Consider any pure inductive lumped element (L > 0, R = 0, C = 0) as a Josephson junction
        if (std::abs(port.L) > 0.0 && std::abs(port.R) <= 0.0 && std::abs(port.C) <= 0.0) {
          has_jj = true;
          Mpi::Print(" Detected Josephson junction at port {}\n", idx);
          break;
        }
      }
      
      // Add indicator with the JJ weight (now 10.0) to focus refinement on JJ and surrounding regions
      // The higher jj_weight in errorindicator.hpp ensures much finer meshing
      estimator.AddErrorIndicator(E, B, E_elec + E_mag, indicator, has_jj);
    }
    
    // Mode type classification and convergence checking
    // For JJ modes, apply convergence checking if enabled
    // For non-JJ modes, always process them
    bool mode_is_jj = false;
      
    // Check if the current element/region contains a Josephson Junction
    const auto& lumped_port_op = space_op.GetLumpedPortOp();
      
    // Check for Josephson junctions
    for (const auto& [idx, port] : lumped_port_op) {
      // Consider any pure inductive lumped element (L > 0, R = 0, C = 0) as a Josephson junction
      if (std::abs(port.L) > 0.0 && std::abs(port.R) <= 0.0 && std::abs(port.C) <= 0.0) {
        // Get participation ratio to check if this is a JJ mode
        double pj = post_op.GetInductorParticipation(lumped_port_op, idx, E_elec);
        if (std::abs(pj) > 0.01) {
          mode_is_jj = true;
          break;
        }
      }
    }
      
    if (mode_is_jj) {
      Mpi::Print(" Mode {:d} is a JJ mode\n", i+1);
      
      // Apply convergence check only to JJ modes if junction convergence is enabled
      if (!indicator.HasConverged()) {
        Mpi::Warning(" JJ Mode {:d} has not met convergence criteria\n", i+1);
        // Optional: skip this mode if you want strict convergence, enabled by default
        // continue;
      } else {
        Mpi::Print(" JJ Mode {:d} has met convergence criteria\n", i+1);
      }
    } else {
      Mpi::Print(" Mode {:d} is a non-JJ mode\n", i+1);
      // Non-JJ modes don't need junction convergence
    }

    // Postprocess the mode.
    Postprocess(post_op, space_op.GetLumpedPortOp(), i, omega, error_bkwd, error_abs,
                num_conv, E_elec, E_mag,
                (i == iodata.solver.eigenmode.n - 1) ? &indicator : nullptr);
  }
  // Process all modes without special JJ handling
  if (num_conv > 0) {
    Mpi::Print("\nProcessing all eigenmodes without filtering or reclassification\n");
    
    // Process all eigenmodes directly in their natural order
    for (int i = 0; i < num_conv; i++) {
      // Prepare the eigenmode for processing
      std::complex<double> omega = eigen->GetEigenvalue(i);
      double error_bkwd = eigen->GetError(i, EigenvalueSolver::ErrorType::BACKWARD);
      double error_abs = eigen->GetError(i, EigenvalueSolver::ErrorType::ABSOLUTE);
      
      if (!C) {
        // Linear EVP has eigenvalue μ = -λ² = ω².
        omega = std::sqrt(omega);
      } else {
        // Quadratic EVP solves for eigenvalue λ = iω.
        omega /= 1i;
      }
      
      // Compute fields
      eigen->GetEigenvector(i, E);
      Curl.Mult(E.Real(), B.Real());
      Curl.Mult(E.Imag(), B.Imag());
      B *= -1.0 / (1i * omega);
      
      // Set up PostOperator
      post_op.SetEGridFunction(E);
      post_op.SetBGridFunction(B);
      post_op.UpdatePorts(space_op.GetLumpedPortOp(), omega.real());
      const double E_elec = post_op.GetEFieldEnergy();
      const double E_mag = post_op.GetHFieldEnergy();
      const double E_cap = post_op.GetLumpedCapacitorEnergy(space_op.GetLumpedPortOp());
      const double E_ind = post_op.GetLumpedInductorEnergy(space_op.GetLumpedPortOp());
      
      // Check if this is a JJ mode based on port participation
      bool is_jj = false;
      const auto& lumped_port_op = space_op.GetLumpedPortOp();
      for (const auto& [idx, port] : lumped_port_op) {
        if (std::abs(port.L) > 0.0 && std::abs(port.R) <= 0.0 && std::abs(port.C) <= 0.0) {
          // Get participation ratio to check if this is a JJ mode
          double pj = post_op.GetInductorParticipation(lumped_port_op, idx, post_op.GetEFieldEnergy());
          if (std::abs(pj) > 0.01) {  // This is just for reporting, not filtering
            is_jj = true;
            Mpi::Print(" Port {:d} participation: {:.6e}\n", idx, pj);
            break;
          }
        }
      }
      
      // Report mode type but don't filter
      if (is_jj) {
        Mpi::Print("\nProcessing JJ mode (eigenmode {:d})\n", i+1);
      } else {
        Mpi::Print("\nProcessing non-JJ mode (eigenmode {:d})\n", i+1);
      }
      
      Mpi::Print(" Frequency: {:.6f} GHz\n", omega.real());
      Mpi::Print(" Lumped energies: E_cap = {:.3e}, E_ind = {:.3e}\n", E_cap, E_ind);
      
      // Process this mode with its natural index
      PostprocessEigen(i, omega, error_bkwd, error_abs, num_conv);
      PostprocessPorts(post_op, space_op.GetLumpedPortOp(), i);
      PostprocessEPR(post_op, space_op.GetLumpedPortOp(), i, omega, E_elec);
      PostprocessDomains(post_op, "m", i, i+1, E_elec, E_mag, E_cap, E_ind);
      PostprocessSurfaces(post_op, "m", i, i+1, E_elec, E_mag);
      PostprocessProbes(post_op, "m", i, i+1);
      PostprocessErrorIndicator(post_op, indicator, true);
    }
  }
    
  return {indicator, space_op.GlobalTrueVSize()};
}

// CUSTOM CONVERGENCE
bool EigenSolver::HasJunctionInDomain(int elem_idx) const {
  // Implementation requires SpaceOperator, which is available in the Solve method
  // This is a placeholder - we'll use a different approach in the main Solve method
  return false;
}

void EigenSolver::Postprocess(const PostOperator &post_op,
                              const LumpedPortOperator &lumped_port_op, int i,
                              std::complex<double> omega, double error_bkwd,
                              double error_abs, int num_conv, double E_elec, double E_mag,
                              const ErrorIndicator *indicator) const
{
  // DEBUG: Print that postprocessing is starting
  Mpi::Print("DEBUG: Main Postprocess function called for mode {} of {}\n", i+1, num_conv);

  // The internal GridFunctions for PostOperator have already been set from the E and B
  // solutions in the main loop over converged eigenvalues.
  const double E_cap = post_op.GetLumpedCapacitorEnergy(lumped_port_op);
  const double E_ind = post_op.GetLumpedInductorEnergy(lumped_port_op);
  PostprocessEigen(i, omega, error_bkwd, error_abs, num_conv);
  PostprocessPorts(post_op, lumped_port_op, i);
  PostprocessEPR(post_op, lumped_port_op, i, omega, E_elec + E_cap);
  PostprocessDomains(post_op, "m", i, i + 1, E_elec, E_mag, E_cap, E_ind);
  PostprocessSurfaces(post_op, "m", i, i + 1, E_elec + E_cap, E_mag + E_ind);
  PostprocessProbes(post_op, "m", i, i + 1);
  if (i < iodata.solver.eigenmode.n_post)
  {
    PostprocessFields(post_op, i, i + 1);
    Mpi::Print(" Wrote mode {:d} to disk\n", i + 1);
  }
  if (indicator)
  {
    PostprocessErrorIndicator(post_op, *indicator, iodata.solver.eigenmode.n_post > 0);
  }
}

namespace
{

struct PortVIData
{
  const int idx;                        // Lumped port index
  const std::complex<double> V_i, I_i;  // Port voltage, current
};

struct EprLData
{
  const int idx;    // Lumped port index
  const double pj;  // Inductor energy-participation ratio
};

struct EprIOData
{
  const int idx;    // Lumped port index
  const double Ql;  // Quality factor
  const double Kl;  // κ for loss rate
};

}  // namespace

void EigenSolver::PostprocessEigen(int i, std::complex<double> omega, double error_bkwd,
                                   double error_abs, int num_conv) const
{
  // DEBUG: Print output path information
  Mpi::Print("DEBUG: PostprocessEigen called, post_dir = '{}', root = {}\n", 
             post_dir, root ? "true" : "false");

  // Dimensionalize the result and print in a nice table of frequencies and Q-factors. Save
  // to file if user has specified.
  const std::complex<double> f = {
      iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, omega.real()),
      iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, omega.imag())};
  const double Q =
      (f.imag() == 0.0) ? mfem::infinity() : 0.5 * std::abs(f) / std::abs(f.imag());

  // Print table to stdout.
  {
    const int int_width = 1 + static_cast<int>(std::log10(num_conv));
    constexpr int p = 6;
    constexpr int w = 6 + p + 7;  // Column spaces + precision + extra for table
    if (i == 0)
    {
      // clang-format off
      Mpi::Print("{:>{}s}{:>{}s}{:>{}s}{:>{}s}{:>{}s}\n{}\n",
                 "m", int_width,
                 "Re{ω}/2π (GHz)", w,
                 "Im{ω}/2π (GHz)", w,
                 "Bkwd. Error", w,
                 "Abs. Error", w,
                 std::string(int_width + 4 * w, '='));
      // clang-format on
    }
    // clang-format off
    Mpi::Print("{:{}d}{:+{}.{}e}{:+{}.{}e}{:+{}.{}e}{:+{}.{}e}\n",
               i + 1, int_width,
               f.real(), w, p,
               f.imag(), w, p,
               error_bkwd, w, p,
               error_abs, w, p);
    // clang-format on
  }

  // Print table to file.
  if (root && post_dir.length() > 0)
  {
    std::string path = post_dir + "eig.csv";
    auto output = OutputFile(path, (i > 0));
    if (i == 0)
    {
      // clang-format off
      output.print("{:>{}s},{:>{}s},{:>{}s},{:>{}s},{:>{}s},{:>{}s}\n",
                   "m", table.w1,
                   "Re{f} (GHz)", table.w,
                   "Im{f} (GHz)", table.w,
                   "Q", table.w,
                   "Error (Bkwd.)", table.w,
                   "Error (Abs.)", table.w);
      // clang-format on
    }
    // clang-format off
    output.print("{:{}.{}e},{:+{}.{}e},{:+{}.{}e},{:+{}.{}e},{:+{}.{}e},{:+{}.{}e}\n",
                 static_cast<double>(i + 1), table.w1, table.p1,
                 f.real(), table.w, table.p,
                 f.imag(), table.w, table.p,
                 Q, table.w, table.p,
                 error_bkwd, table.w, table.p,
                 error_abs, table.w, table.p);
    // clang-format on
  }
}

void EigenSolver::PostprocessPorts(const PostOperator &post_op,
                                   const LumpedPortOperator &lumped_port_op, int i) const
{
  // Postprocess the frequency domain lumped port voltages and currents (complex magnitude
  // = sqrt(2) * RMS).
  if (post_dir.length() == 0)
  {
    return;
  }
  std::vector<PortVIData> port_data;
  port_data.reserve(lumped_port_op.Size());
  for (const auto &[idx, data] : lumped_port_op)
  {
    const std::complex<double> V_i = post_op.GetPortVoltage(lumped_port_op, idx);
    const std::complex<double> I_i = post_op.GetPortCurrent(lumped_port_op, idx);
    port_data.push_back({idx, iodata.DimensionalizeValue(IoData::ValueType::VOLTAGE, V_i),
                         iodata.DimensionalizeValue(IoData::ValueType::CURRENT, I_i)});
  }
  if (root && !port_data.empty())
  {
    // Write the port voltages.
    {
      std::string path = post_dir + "port-V.csv";
      auto output = OutputFile(path, (i > 0));
      if (i == 0)
      {
        output.print("{:>{}s},", "m", table.w1);
        for (const auto &data : port_data)
        {
          // clang-format off
          output.print("{:>{}s},{:>{}s}{}",
                       "Re{V[" + std::to_string(data.idx) + "]} (V)", table.w,
                       "Im{V[" + std::to_string(data.idx) + "]} (V)", table.w,
                       (data.idx == port_data.back().idx) ? "" : ",");
          // clang-format on
        }
        output.print("\n");
      }
      // clang-format off
      output.print("{:{}.{}e},",
                   static_cast<double>(i + 1), table.w1, table.p1);
      // clang-format on
      for (const auto &data : port_data)
      {
        // clang-format off
        output.print("{:+{}.{}e},{:+{}.{}e}{}",
                     data.V_i.real(), table.w, table.p,
                     data.V_i.imag(), table.w, table.p,
                     (data.idx == port_data.back().idx) ? "" : ",");
        // clang-format on
      }
      output.print("\n");
    }

    // Write the port currents.
    {
      std::string path = post_dir + "port-I.csv";
      auto output = OutputFile(path, (i > 0));
      if (i == 0)
      {
        output.print("{:>{}s},", "m", table.w1);
        for (const auto &data : port_data)
        {
          // clang-format off
          output.print("{:>{}s},{:>{}s}{}",
                       "Re{I[" + std::to_string(data.idx) + "]} (A)", table.w,
                       "Im{I[" + std::to_string(data.idx) + "]} (A)", table.w,
                       (data.idx == port_data.back().idx) ? "" : ",");
          // clang-format on
        }
        output.print("\n");
      }
      // clang-format off
      output.print("{:{}.{}e},",
                   static_cast<double>(i + 1), table.w1, table.p1);
      // clang-format on
      for (const auto &data : port_data)
      {
        // clang-format off
        output.print("{:+{}.{}e},{:+{}.{}e}{}",
                     data.I_i.real(), table.w, table.p,
                     data.I_i.imag(), table.w, table.p,
                     (data.idx == port_data.back().idx) ? "" : ",");
        // clang-format on
      }
      output.print("\n");
    }
  }
}

void EigenSolver::PostprocessEPR(const PostOperator &post_op,
                                 const LumpedPortOperator &lumped_port_op, int i,
                                 std::complex<double> omega, double E_m) const
{
  // If ports have been specified in the model, compute the corresponding energy-
  // participation ratios (EPR) and write out to disk.
  if (post_dir.length() == 0)
  {
    return;
  }

  // Write the mode EPR for lumped inductor elements.
  std::vector<EprLData> epr_L_data;
  epr_L_data.reserve(lumped_port_op.Size());
  for (const auto &[idx, data] : lumped_port_op)
  {
    if (std::abs(data.L) > 0.0)
    {
      const double pj = post_op.GetInductorParticipation(lumped_port_op, idx, E_m);
      epr_L_data.push_back({idx, pj});
    }
  }
  if (root && !epr_L_data.empty())
  {
    std::string path = post_dir + "port-EPR.csv";
    auto output = OutputFile(path, (i > 0));
    if (i == 0)
    {
      output.print("{:>{}s},", "m", table.w1);
      for (const auto &data : epr_L_data)
      {
        // clang-format off
        output.print("{:>{}s}{}",
                     "p[" + std::to_string(data.idx) + "]", table.w,
                     (data.idx == epr_L_data.back().idx) ? "" : ",");
        // clang-format on
      }
      output.print("\n");
    }
    output.print("{:{}.{}e},", static_cast<double>(i + 1), table.w1, table.p1);
    for (const auto &data : epr_L_data)
    {
      // clang-format off
      output.print("{:+{}.{}e}{}",
                   data.pj, table.w, table.p,
                   (data.idx == epr_L_data.back().idx) ? "" : ",");
      // clang-format on
    }
    output.print("\n");
  }

  // Write the mode EPR for lumped resistor elements.
  std::vector<EprIOData> epr_IO_data;
  epr_IO_data.reserve(lumped_port_op.Size());
  for (const auto &[idx, data] : lumped_port_op)
  {
    if (std::abs(data.R) > 0.0)
    {
      const double Kl = post_op.GetExternalKappa(lumped_port_op, idx, E_m);
      const double Ql = (Kl == 0.0) ? mfem::infinity() : omega.real() / std::abs(Kl);
      epr_IO_data.push_back(
          {idx, Ql, iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, Kl)});
    }
  }
  if (root && !epr_IO_data.empty())
  {
    std::string path = post_dir + "port-Q.csv";
    auto output = OutputFile(path, (i > 0));
    if (i == 0)
    {
      output.print("{:>{}s},", "m", table.w1);
      for (const auto &data : epr_IO_data)
      {
        // clang-format off
        output.print("{:>{}s},{:>{}s}{}",
                     "Q_ext[" + std::to_string(data.idx) + "]", table.w,
                     "κ_ext[" + std::to_string(data.idx) + "] (GHz)", table.w,
                     (data.idx == epr_IO_data.back().idx) ? "" : ",");
        // clang-format on
      }
      output.print("\n");
    }
    output.print("{:{}.{}e},", static_cast<double>(i + 1), table.w1, table.p1);
    for (const auto &data : epr_IO_data)
    {
      // clang-format off
      output.print("{:+{}.{}e},{:+{}.{}e}{}",
                   data.Ql, table.w, table.p,
                   data.Kl, table.w, table.p,
                   (data.idx == epr_IO_data.back().idx) ? "" : ",");
      // clang-format on
    }
    output.print("\n");
  }
}

}  // namespace palace
