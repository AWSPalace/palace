// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "magnetostaticsolver.hpp"

#include <mfem.hpp>
#include "fem/errorindicator.hpp"
#include "fem/mesh.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "models/curlcurloperator.hpp"
#include "models/postoperator.hpp"
#include "models/surfacecurrentoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"

namespace palace
{

std::pair<ErrorIndicator, long long int>
MagnetostaticSolver::Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const
{
  // Construct the system matrix defining the linear operator. Dirichlet boundaries are
  // handled eliminating the rows and columns of the system matrix for the corresponding
  // dofs.
  BlockTimer bt0(Timer::CONSTRUCT);
  CurlCurlOperator curlcurlop(iodata, mesh);
  auto K = curlcurlop.GetStiffnessMatrix();
  const auto &Curl = curlcurlop.GetCurlMatrix();
  SaveMetadata(curlcurlop.GetNDSpaces());

  // Set up the linear solver.
  KspSolver ksp(iodata, curlcurlop.GetNDSpaces(), &curlcurlop.GetH1Spaces());
  ksp.SetOperators(*K, *K);

  // Terminal indices are the set of boundaries over which to compute the inductance matrix.
  PostOperator postop(iodata, curlcurlop, "magnetostatic");
  int nstep = static_cast<int>(curlcurlop.GetSurfaceCurrentOp().Size());
  MFEM_VERIFY(nstep > 0,
              "No surface current boundaries specified for magnetostatic simulation!");

  // Source term and solution vector storage.
  Vector RHS(Curl.Width()), B(Curl.Height());
  std::vector<Vector> A(nstep);
  std::vector<double> I_inc(nstep);

  // Initialize structures for storing and reducing the results of error estimation.
  CurlFluxErrorEstimator estimator(
      curlcurlop.GetMaterialOp(), curlcurlop.GetRTSpace(), curlcurlop.GetNDSpaces(),
      iodata.solver.linear.estimator_tol, iodata.solver.linear.estimator_max_it, 0,
      iodata.solver.linear.estimator_mg);
  ErrorIndicator indicator;

  // Main loop over current source boundaries.
  Mpi::Print("\nComputing magnetostatic fields for {:d} source boundar{}\n", nstep,
             (nstep > 1) ? "ies" : "y");
  int step = 0;
  auto t0 = Timer::Now();
  for (const auto &[idx, data] : curlcurlop.GetSurfaceCurrentOp())
  {
    Mpi::Print("\nIt {:d}/{:d}: Index = {:d} (elapsed time = {:.2e} s)\n", step + 1, nstep,
               idx, Timer::Duration(Timer::Now() - t0).count());

    // Form and solve the linear system for a prescribed current on the specified source.
    Mpi::Print("\n");
    A[step].SetSize(RHS.Size());
    A[step].UseDevice(true);
    A[step] = 0.0;
    curlcurlop.GetExcitationVector(idx, RHS);
    ksp.Mult(RHS, A[step]);

    // Compute B = ∇ x A on the true dofs, and set the internal GridFunctions in
    // PostOperator for all postprocessing operations.
    BlockTimer bt2(Timer::POSTPRO);
    Curl.Mult(A[step], B);
    postop.SetAGridFunction(A[step]);
    postop.SetBGridFunction(B);
    const double E_mag = postop.GetHFieldEnergy();
    Mpi::Print(" Sol. ||A|| = {:.6e} (||RHS|| = {:.6e})\n",
               linalg::Norml2(curlcurlop.GetComm(), A[step]),
               linalg::Norml2(curlcurlop.GetComm(), RHS));
    {
      const double J = iodata.DimensionalizeValue(IoData::ValueType::ENERGY, 1.0);
      Mpi::Print(" Field energy H = {:.3e} J\n", E_mag * J);
    }
    I_inc[step] = data.GetExcitationCurrent();

    // Calculate and record the error indicators.
    Mpi::Print(" Updating solution error estimates\n");
    estimator.AddErrorIndicator(B, E_mag, indicator);

    // Postprocess field solutions and optionally write solution to disk.
    Postprocess(postop, step, idx, I_inc[step], E_mag,
                (step == nstep - 1) ? &indicator : nullptr);

    // Next source.
    step++;
  }

  // Postprocess the inductance matrix from the computed field solutions.
  BlockTimer bt1(Timer::POSTPRO);
  SaveMetadata(ksp);
  PostprocessTerminals(postop, curlcurlop.GetSurfaceCurrentOp(), A, I_inc);
  return {indicator, curlcurlop.GlobalTrueVSize()};
}

void MagnetostaticSolver::Postprocess(const PostOperator &postop, int step, int idx,
                                      double I_inc, double E_mag,
                                      const ErrorIndicator *indicator) const
{
  // The internal GridFunctions for PostOperator have already been set from the A solution
  // in the main loop.
  PostprocessDomains(postop, "i", step, idx, 0.0, E_mag, 0.0, 0.0);
  PostprocessSurfaces(postop, "i", step, idx, 0.0, E_mag);
  PostprocessProbes(postop, "i", step, idx);
  if (step < iodata.solver.magnetostatic.n_post)
  {
    PostprocessFields(postop, step, idx);
    Mpi::Print(" Wrote fields to disk for source {:d}\n", idx);
  }
  if (indicator)
  {
    PostprocessErrorIndicator(postop, *indicator, iodata.solver.magnetostatic.n_post > 0);
  }
}

void MagnetostaticSolver::PostprocessTerminals(PostOperator &postop,
                                               const SurfaceCurrentOperator &surf_j_op,
                                               const std::vector<Vector> &A,
                                               const std::vector<double> &I_inc) const
{
  // Postprocess the Maxwell inductance matrix. See p. 97 of the COMSOL AC/DC Module manual
  // for the associated formulas based on the magnetic field energy based on a current
  // excitation for each port. Alternatively, we could compute the resulting loop fluxes to
  // get M directly as:
  //                         Φ_i = ∫ B ⋅ n_j dS
  // and M_ij = Φ_i/I_j. The energy formulation avoids having to locally integrate B =
  // ∇ x A.
  mfem::DenseMatrix M(A.size()), Mm(A.size());
  for (int i = 0; i < M.Height(); i++)
  {
    // Diagonal: Mᵢᵢ = 2 Uₘ(Aᵢ) / Iᵢ² = (Aᵢᵀ K Aᵢ) / Iᵢ²
    auto &A_gf = postop.GetAGridFunction().Real();
    auto &H_gf = postop.GetDomainPostOp().H;
    A_gf.SetFromTrueDofs(A[i]);
    postop.GetDomainPostOp().M_mag->Mult(A_gf, H_gf);
    M(i, i) = Mm(i, i) =
        linalg::Dot<Vector>(postop.GetComm(), A_gf, H_gf) / (I_inc[i] * I_inc[i]);

    // Off-diagonals: Mᵢⱼ = Uₘ(Aᵢ + Aⱼ) / (Iᵢ Iⱼ) - 1/2 (Iᵢ/Iⱼ Mᵢᵢ + Iⱼ/Iᵢ Mⱼⱼ)
    //                    = (Aⱼᵀ K Aᵢ) / (Iᵢ Iⱼ)
    for (int j = i + 1; j < M.Width(); j++)
    {
      A_gf.SetFromTrueDofs(A[j]);
      M(i, j) = linalg::Dot<Vector>(postop.GetComm(), A_gf, H_gf) / (I_inc[i] * I_inc[j]);
      Mm(i, j) = -M(i, j);
      Mm(i, i) -= Mm(i, j);
    }

    // Copy lower triangle from already computed upper triangle.
    for (int j = 0; j < i; j++)
    {
      M(i, j) = M(j, i);
      Mm(i, j) = Mm(j, i);
      Mm(i, i) -= Mm(i, j);
    }
  }
  mfem::DenseMatrix Minv(M);
  Minv.Invert();  // In-place, uses LAPACK (when available) and should be cheap

  // Only root writes to disk (every process has full matrices).
  if (!root || post_dir.length() == 0)
  {
    return;
  }

  // Write inductance matrix data.
  auto PrintMatrix = [&surf_j_op, this](const std::string &file, const std::string &name,
                                        const std::string &unit,
                                        const mfem::DenseMatrix &mat, double scale)
  {
    std::string path = post_dir + file;
    auto output = OutputFile(path, false);
    output.print("{:>{}s},", "i", table.w1);
    for (const auto &[idx2, data2] : surf_j_op)
    {
      // clang-format off
      output.print("{:>{}s}{}",
                   name + "[i][" + std::to_string(idx2) + "] " + unit, table.w,
                   (idx2 == surf_j_op.rbegin()->first) ? "" : ",");
      // clang-format on
    }
    output.print("\n");
    int i = 0;
    for (const auto &[idx, data] : surf_j_op)
    {
      int j = 0;
      output.print("{:{}.{}e},", static_cast<double>(idx), table.w1, table.p1);
      for (const auto &[idx2, data2] : surf_j_op)
      {
        // clang-format off
        output.print("{:+{}.{}e}{}",
                     mat(i, j) * scale, table.w, table.p,
                     (idx2 == surf_j_op.rbegin()->first) ? "" : ",");
        // clang-format on
        j++;
      }
      output.print("\n");
      i++;
    }
  };
  const double H = iodata.DimensionalizeValue(IoData::ValueType::INDUCTANCE, 1.0);
  PrintMatrix("terminal-M.csv", "M", "(H)", M, H);
  PrintMatrix("terminal-Minv.csv", "M⁻¹", "(1/H)", Minv, 1.0 / H);
  PrintMatrix("terminal-Mm.csv", "M_m", "(H)", Mm, H);

  // Also write out a file with source current excitations.
  {
    std::string path = post_dir + "terminal-I.csv";
    auto output = OutputFile(path, false);
    // clang-format off
    output.print("{:>{}s},{:>{}s}\n",
                 "i", table.w1,
                 "I_inc[i] (A)", table.w);
    // clang-format on
    int i = 0;
    for (const auto &[idx, data] : surf_j_op)
    {
      // clang-format off
      output.print("{:{}.{}e},{:+{}.{}e}\n",
                   static_cast<double>(idx), table.w1, table.p1,
                   iodata.DimensionalizeValue(IoData::ValueType::CURRENT, I_inc[i]),
                   table.w, table.p);
      // clang-format on
      i++;
    }
  }
}

}  // namespace palace
