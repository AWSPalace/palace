// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "errorestimator.hpp"
#include "fem/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/errorindicators.hpp"
#include "utils/iodata.hpp"
#include "utils/mfemcoefficients.hpp"
#include "utils/mfemintegrators.hpp"
#include "utils/multigrid.hpp"

namespace palace
{

namespace
{

// Given a grid function defining a vector solution, compute the error relative
// to a vector coefficient. The default quadrature rule exactly integrates 2p +
// q polynomials, but the quadrature order can be increased or decreased via
// quad_order_increment.
mfem::Vector ComputeElementLpErrors(const mfem::ParGridFunction &sol, double p,
                                    mfem::VectorCoefficient &exsol,
                                    int quad_order_increment = 0);

// Given a grid function defining a vector/scalar solution, compute the Lp norm of the
// solution. The default quadrature rule exactly integrates 2p +
// q polynomials, but the quadrature order can be increased or decreased via
// quad_order_increment for non L2 norms.
double ComputeVectorLpNorm(const mfem::ParGridFunction &sol, double p = 2,
                           int quad_order_increment = 0);
double ComputeScalarLpNorm(const mfem::ParGridFunction &sol, double p = 2,
                           int quad_order_increment = 0);

}  // namespace

using namespace utils;

CurlFluxErrorEstimator::CurlFluxErrorEstimator(
    const IoData &iodata, const MaterialOperator &mat_op,
    const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh,
    mfem::ParFiniteElementSpace &fes)
  : mat_op(mat_op), fes(fes), flux_fec(iodata.solver.order - 1, mesh.back()->Dimension()),
    flux_fes(mesh.back().get(), &flux_fec),
    smooth_flux_fecs(ConstructFECollections<mfem::ND_FECollection>(
        iodata.solver.linear.mat_gmg, false, iodata.solver.order,
        mesh.back()->Dimension())),
    smooth_flux_fes(utils::ConstructFiniteElementSpaceHierarchy<mfem::ND_FECollection>(
        mesh, smooth_flux_fecs)),
    projector(smooth_flux_fes, iodata.solver.linear.tol / 10, 200, iodata.problem.verbose)
{
}

CurlFluxErrorEstimator::CurlFluxErrorEstimator(const IoData &iodata,
                                               const MaterialOperator &mat_op,
                                               std::unique_ptr<mfem::ParMesh> &mesh,
                                               mfem::ParFiniteElementSpace &fes)
  : mat_op(mat_op), fes(fes), flux_fec(iodata.solver.order - 1, mesh->Dimension()),
    flux_fes(mesh.get(), &flux_fec),
    smooth_flux_fecs(ConstructFECollections<mfem::ND_FECollection>(
        iodata.solver.linear.mat_gmg, false, iodata.solver.order, mesh->Dimension())),
    smooth_flux_fes(utils::ConstructFiniteElementSpaceHierarchy<mfem::ND_FECollection>(
        mesh, smooth_flux_fecs)),
    projector(smooth_flux_fes, iodata.solver.linear.tol / 10, 200, iodata.problem.verbose)
{
}

mfem::Vector CurlFluxErrorEstimator::operator()(const petsc::PetscParVector &v,
                                                bool use_mfem) const
{
  const auto cv = v.GetToVectors();
  mfem::ParComplexGridFunction field(&fes);
  field.real().SetFromTrueDofs(cv.real);
  field.imag().SetFromTrueDofs(cv.imag);

  const int nelem = smooth_flux_fes.GetFinestFESpace().GetNE();
  mfem::Vector real_error(nelem), imag_error(nelem);

  constexpr int normp = 2;  // 2 norm ensures no under integration.
  if (use_mfem)
  {
    // This is used for comparison purposes only.
    // TODO: Delete this once the alternative direct construction is ready.
    // NB: There are bugs within L2ZZErrorEstimator presently, that ND simplices
    // with p > 1 will not work, and also a flux fec of RT elements will also
    // give incorrect answers.
    const int order = fes.GetElementOrder(0);  // Assumes no mixed p.
    auto *const pmesh = fes.GetParMesh();

    MaterialPropertyCoefficient<MaterialPropertyType::INV_PERMEABILITY> coef(mat_op);
    mfem::CurlCurlIntegrator flux_integrator(coef);  // L2ZZErrorEstimator ignores coef

    static_cast<void>(mfem::L2ZZErrorEstimator(flux_integrator, field.real(),
                                               smooth_flux_fes.GetFinestFESpace(), flux_fes,
                                               real_error, normp));

    static_cast<void>(mfem::L2ZZErrorEstimator(flux_integrator, field.imag(),
                                               smooth_flux_fes.GetFinestFESpace(), flux_fes,
                                               imag_error, normp));
  }
  else
  {
    // This lambda computes the flux as achieved within mfem. Use this to bench against.
    auto mfem_flux_func = [this, &field]()
    {
      mfem::ParComplexGridFunction flux_func(&flux_fes);
      flux_func.real() = 0.0;
      flux_func.imag() = 0.0;

      mfem::Array<int> xdofs, fdofs;
      mfem::Vector el_x, el_f;
      mfem::CurlCurlIntegrator curl;

      for (int i = 0; i < fes.GetNE(); ++i)
      {
        const auto *const xdoftrans = fes.GetElementVDofs(i, xdofs);
        field.real().GetSubVector(xdofs, el_x);
        if (xdoftrans)
        {
          xdoftrans->InvTransformPrimal(el_x);
        }

        auto *T = field.real().ParFESpace()->GetElementTransformation(i);
        curl.ComputeElementFlux(*field.real().ParFESpace()->GetFE(i), *T, el_x,
                                *flux_fes.GetFE(i), el_f, false);

        const auto *const fdoftrans = flux_fes.GetElementVDofs(i, fdofs);
        if (fdoftrans)
        {
          fdoftrans->TransformPrimal(el_f);
        }

        flux_func.real().SetSubVector(fdofs, el_f);

        field.imag().GetSubVector(xdofs, el_x);
        if (xdoftrans)
        {
          xdoftrans->InvTransformPrimal(el_x);
        }
        curl.ComputeElementFlux(*field.imag().ParFESpace()->GetFE(i), *T, el_x,
                                *flux_fes.GetFE(i), el_f, false);
        if (fdoftrans)
        {
          fdoftrans->TransformPrimal(el_f);
        }
        flux_func.imag().SetSubVector(fdofs, el_f);
      }

      return flux_func;
    };

    // This lambda computes the curl flux function by forming a projector. In
    // theory same as the above, but still doesn't allow for coefficients. This
    // is another benching device, helpful for debugging.
    auto projector_flux_func = [this, &cv, &field]()
    {
      // Interpolate the weighted curl of the field in fes onto the discrete
      // flux space.
      mfem::ParDiscreteLinearOperator curl(&fes, &flux_fes);  // (domain, range)

      curl.AddDomainInterpolator(new mfem::CurlInterpolator);
      curl.Assemble();
      curl.Finalize();
      auto hyprecurlop = std::unique_ptr<mfem::HypreParMatrix>(curl.ParallelAssemble());

      ComplexVector flux(flux_fes.GetTrueVSize());

      hyprecurlop->Mult(cv.real, flux.real);
      hyprecurlop->Mult(cv.imag, flux.imag);

      mfem::ParComplexGridFunction flux_func(&flux_fes);
      flux_func.real().SetFromTrueDofs(flux.real);
      flux_func.imag().SetFromTrueDofs(flux.imag);

      return flux_func;
    };

    // Given a flux function built elsewhere, construct the linear operator RHS.
    auto flux_rhs_from_func = [&v](auto &fes, mfem::ParComplexGridFunction &flux_func)
    {
      const auto ndof = fes.GetTrueVSize();

      ComplexVector flux(ndof);
      {
        mfem::VectorGridFunctionCoefficient f_real(&flux_func.real());
        mfem::ParLinearForm rhs(&fes);
        rhs.AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(f_real));
        rhs.UseFastAssembly(true);
        rhs.Assemble();
        rhs.ParallelAssemble(flux.real);
      }

      {
        mfem::VectorGridFunctionCoefficient f_imag(&flux_func.imag());
        mfem::ParLinearForm rhs(&fes);
        rhs.AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(f_imag));
        rhs.UseFastAssembly(true);
        rhs.Assemble();
        rhs.ParallelAssemble(flux.imag);
      }

      return flux;
    };

    // Compare the two flux functions, should be identical.
    auto comp = [this](const auto &mflux, const auto &pflux)
    {
      bool match = true;

      mfem::Array<int> flux_dofs;
      mfem::Vector mflux_val, pflux_val;
      for (int e = 0; e < flux_fes.GetNE(); e++)
      {
        flux_fes.GetElementVDofs(e, flux_dofs);

        pflux.GetSubVector(flux_dofs, pflux_val);
        mflux.GetSubVector(flux_dofs, mflux_val);

        MFEM_ASSERT(mflux_val.Size() == pflux_val.Size(), "Must match size.");
        constexpr double tol = 1e-6;
        for (int i = 0; i < mflux_val.Size(); ++i)
        {
          auto diff = std::abs(mflux_val(i) - pflux_val(i));
          if (diff > tol)
          {
            std::cout << "Mismatch on e " << e << " i " << i << ": " << mflux_val(i) << " "
                      << pflux_val(i) << '\n';
            match = false;
          }
        }
      }

      return match;
    };

    // MFEM_ASSERT(comp(mfem_flux_func().real(), projector_flux_func().real()),
    //             "Mismatch between projector and L2ZZ construction real values");
    // MFEM_ASSERT(comp(mfem_flux_func().imag(), projector_flux_func().imag()),
    //             "Mismatch between projector and L2ZZ construction imag values");

    // Coefficients for computing the discontinuous flux., i.e. (W, μ⁻¹∇ × V).
    // The code from here down will ultimately be the way to calculate the flux.
    CurlFluxCoefficient real_coef(field.real(), mat_op), imag_coef(field.imag(), mat_op);
    auto rhs_from_coef = [](mfem::ParFiniteElementSpace &smooth_flux_fes, auto &coef)
    {
      mfem::Vector RHS(smooth_flux_fes.GetTrueVSize());

      mfem::ParLinearForm rhs(&smooth_flux_fes);
      rhs.AddDomainIntegrator(new VectorFEDomainLFIntegrator(coef));
      rhs.UseFastAssembly(true);
      rhs.Assemble();
      rhs.ParallelAssemble(RHS);

      return RHS;
    };

    // Switching between these two gives identical.
    // const auto flux = flux_rhs_from_func(smooth_flux_fes.GetFinestFESpace(), projector_flux_func());
    const auto flux =
        ComplexVector(rhs_from_coef(smooth_flux_fes.GetFinestFESpace(), real_coef),
                      rhs_from_coef(smooth_flux_fes.GetFinestFESpace(), imag_coef));

    const auto pflux = petsc::PetscParVector(v.GetComm(), flux);

    // Given the RHS vector of non-smooth flux, construct a flux projector and
    // perform mass matrix inversion in the appropriate space, giving f = M⁻¹ f̂.
    auto build_smooth_flux = [this](const petsc::PetscParVector &flux)
    {
      // Use a copy construction to match appropriate size.
      petsc::PetscParVector smooth_flux(flux);
      projector.Mult(flux, smooth_flux);
      return smooth_flux;
    };
    auto smooth_flux = build_smooth_flux(pflux);

    // Given a complex solution represented with a PetscParVector, build a
    // ComplexGridFunction for evaluation.
    auto build_func = [](const petsc::PetscParVector &f, mfem::ParFiniteElementSpace &fes)
    {
      mfem::ParComplexGridFunction flux(&fes);
      const auto fi = f.GetToVectors();
      flux.real().SetFromTrueDofs(fi.real);
      flux.imag().SetFromTrueDofs(fi.imag);
      flux.real().ExchangeFaceNbrData();
      flux.imag().ExchangeFaceNbrData();
      return flux;
    };

    auto smooth_flux_func = build_func(smooth_flux, smooth_flux_fes.GetFinestFESpace());

    real_error = ComputeElementLpErrors(smooth_flux_func.real(), normp, real_coef);
    imag_error = ComputeElementLpErrors(smooth_flux_func.imag(), normp, imag_coef);
  }

  // Compute the magnitude of the complex valued error.
  auto magnitude = [](const auto &r, const auto &i) { return std::sqrt(r * r + i * i); };
  mfem::Vector estimates(real_error.Size());
  std::transform(real_error.begin(), real_error.end(), imag_error.begin(),
                 estimates.begin(), magnitude);

  // Normalize the error by the solution L2 norm.
  const auto normalization = std::sqrt(std::pow(ComputeVectorLpNorm(field.real(), 2), 2.0) +
                                       std::pow(ComputeVectorLpNorm(field.imag(), 2), 2.0));

  std::for_each(estimates.begin(), estimates.end(),
                [&normalization](auto &x) { x /= normalization; });

  return estimates;
}

GradFluxErrorEstimator::GradFluxErrorEstimator(
    const IoData &iodata, const MaterialOperator &mat_op,
    const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh,
    mfem::ParFiniteElementSpace &fes)
  : mat_op(mat_op), fes(fes), flux_fec(iodata.solver.order - 1, mesh.back()->Dimension()),
    flux_fes(mesh.back().get(), &flux_fec, mesh.back()->Dimension()),
    smooth_flux_fecs(ConstructFECollections<mfem::H1_FECollection>(
        iodata.solver.linear.mat_gmg, false, iodata.solver.order,
        mesh.back()->Dimension())),
    smooth_flux_fes(utils::ConstructFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
        mesh, smooth_flux_fecs, mesh.back()->Dimension())),
    projector(smooth_flux_fes, iodata.solver.linear.tol / 10, 200, iodata.problem.verbose)
{
}

GradFluxErrorEstimator::GradFluxErrorEstimator(const IoData &iodata,
                                               const MaterialOperator &mat_op,
                                               std::unique_ptr<mfem::ParMesh> &mesh,
                                               mfem::ParFiniteElementSpace &fes)
  : mat_op(mat_op), fes(fes), flux_fec(iodata.solver.order - 1, mesh->Dimension()),
    flux_fes(mesh.get(), &flux_fec, mesh->Dimension()),
    smooth_flux_fecs(ConstructFECollections<mfem::H1_FECollection>(
        iodata.solver.linear.mat_gmg, false, iodata.solver.order, mesh->Dimension())),
    smooth_flux_fes(utils::ConstructFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
        mesh, smooth_flux_fecs, mesh->Dimension())),
    projector(smooth_flux_fes, iodata.solver.linear.tol / 10, 200, iodata.problem.verbose)
{
}

mfem::Vector GradFluxErrorEstimator::operator()(const mfem::Vector &v, bool use_mfem) const
{
  mfem::ParGridFunction field(&fes);
  field.SetFromTrueDofs(v);

  const int nelem = smooth_flux_fes.GetFinestFESpace().GetNE();

  constexpr int normp = 2;  // 2 norm ensures no under integration.
  if (use_mfem)
  {
    mfem::Vector error(nelem);
    // This is used for comparison purposes only.
    // TODO: Delete this once the alternative direct construction is ready.
    const int order = fes.GetElementOrder(0);  // Assumes no mixed p.
    auto *const pmesh = fes.GetParMesh();

    MaterialPropertyCoefficient<MaterialPropertyType::PERMITTIVITY_ABS> eps(mat_op);
    mfem::DiffusionIntegrator flux_integrator(eps);

    static_cast<void>(mfem::L2ZZErrorEstimator(flux_integrator, field,
                                               smooth_flux_fes.GetFinestFESpace(), flux_fes,
                                               error, normp));

    // Normalize the error by the solution L2 norm, ensures reductions are well scaled.
    const auto normalization = ComputeScalarLpNorm(field, 2);
    std::for_each(error.begin(), error.end(),
                  [&normalization](auto &x) { x /= normalization; });
    return error;
  }
  else
  {
    // This lambda computes the flux as achieved within mfem. Use this to bench against.
    auto mfem_flux_func = [this, &field]()
    {
      mfem::ParGridFunction flux_func(&flux_fes);
      flux_func = 0.0;

      mfem::Array<int> xdofs, fdofs;
      mfem::Vector el_x, el_f;
      MaterialPropertyCoefficient<MaterialPropertyType::PERMITTIVITY_ABS> eps(mat_op);
      mfem::DiffusionIntegrator grad(eps);

      for (int i = 0; i < fes.GetNE(); ++i)
      {
        const auto *const xdoftrans = fes.GetElementVDofs(i, xdofs);
        field.GetSubVector(xdofs, el_x);
        if (xdoftrans)
        {
          xdoftrans->InvTransformPrimal(el_x);
        }

        auto *T = field.ParFESpace()->GetElementTransformation(i);
        grad.ComputeElementFlux(*field.ParFESpace()->GetFE(i), *T, el_x,
        *flux_fes.GetFE(i),
                                el_f, false);

        const auto *const fdoftrans = flux_fes.GetElementVDofs(i, fdofs);
        if (fdoftrans)
        {
          fdoftrans->TransformPrimal(el_f);
        }

        flux_func.SetSubVector(fdofs, el_f);
      }

      return flux_func;
    };

    // This lambda computes the grad flux function by forming a projector. In
    // theory same as the above, but still doesn't allow for coefficients. This
    // is another benching device, helpful for debugging.
    auto projector_flux_func = [this, &v, &field]()
    {
      // Interpolate the weighted grad of the field in fes onto the discrete
      // flux space.
      mfem::ParDiscreteLinearOperator grad(&fes, &flux_fes);  // (domain, range)

      grad.AddDomainInterpolator(new mfem::GradientInterpolator);
      grad.Assemble();
      grad.Finalize();
      auto hypregradop = std::unique_ptr<mfem::HypreParMatrix>(grad.ParallelAssemble());

      mfem::Vector flux(flux_fes.GetTrueVSize());

      hypregradop->Mult(v, flux);

      mfem::ParGridFunction flux_func(&flux_fes);
      flux_func.SetFromTrueDofs(flux);

      return flux_func;
    };

    // Given a flux function built elsewhere, construct the linear operator RHS.
    auto flux_rhs_from_func = [&v](auto &fes, mfem::ParGridFunction &flux_func)
    {
      mfem::Vector flux(fes.GetTrueVSize());

      mfem::VectorGridFunctionCoefficient f(&flux_func);
      mfem::ParLinearForm rhs(&fes);
      rhs.AddDomainIntegrator(new mfem::VectorDomainLFIntegrator(f));
      rhs.UseFastAssembly(true);
      rhs.Assemble();
      rhs.ParallelAssemble(flux);

      return flux;
    };

    // Compare the two flux functions, should be identical.
    auto comp = [this](const auto &mflux, const auto &pflux)
    {
      bool match = true;

      mfem::Array<int> flux_dofs;
      mfem::Vector mflux_val, pflux_val;
      for (int e = 0; e < flux_fes.GetNE(); e++)
      {
        flux_fes.GetElementVDofs(e, flux_dofs);

        pflux.GetSubVector(flux_dofs, pflux_val);
        mflux.GetSubVector(flux_dofs, mflux_val);

        MFEM_ASSERT(mflux_val.Size() == pflux_val.Size(), "Must match size.");
        constexpr double tol = 1e-6;
        for (int i = 0; i < mflux_val.Size(); ++i)
        {
          auto diff = std::abs(mflux_val(i) - pflux_val(i));
          if (diff > tol)
          {
            std::cout << "Mismatch on e " << e << " i " << i << ": " << mflux_val(i) << "
            "
                      << pflux_val(i) << '\n';
            match = false;
          }
        }
      }

      return match;
    };

    // MFEM_ASSERT(comp(mfem_flux_func(), projector_flux_func()),
    //             "Mismatch between projector and L2ZZ construction values");

    // Coefficients for computing the discontinuous flux., i.e. (V, ϵ ∇ ϕ).
    // The code from here down will ultimately be the way to calculate the flux.
    GradFluxCoefficient coef(field, mat_op);
    auto rhs_from_coef = [](mfem::ParFiniteElementSpace &smooth_flux_fes, auto &coef)
    {
      mfem::Vector RHS(smooth_flux_fes.GetTrueVSize());

      mfem::ParLinearForm rhs(&smooth_flux_fes);
      rhs.AddDomainIntegrator(new mfem::VectorDomainLFIntegrator(coef));
      rhs.UseFastAssembly(true);
      rhs.Assemble();
      rhs.ParallelAssemble(RHS);

      return RHS;
    };

    // Switching between these two gives identical.
    // auto flux_func = mfem_flux_func();
    // const auto flux = flux_rhs_from_func(smooth_flux_fes.GetFinestFESpace(), flux_func);
    const auto flux = rhs_from_coef(smooth_flux_fes.GetFinestFESpace(), coef);

    // Given the RHS vector of non-smooth flux, construct a flux projector and
    // perform mass matrix inversion in the appropriate space, giving f = M⁻¹ f̂.
    auto build_smooth_flux = [this](const mfem::Vector &flux)
    {
      // Use a copy construction to match appropriate size.
      mfem::Vector smooth_flux(flux);
      projector.Mult(flux, smooth_flux);
      return smooth_flux;
    };
    auto smooth_flux = build_smooth_flux(flux);

    // Given a solution represented with a Vector, build a GridFunction for evaluation.
    auto build_func = [](const mfem::Vector &f, mfem::ParFiniteElementSpace &fes)
    {
      mfem::ParGridFunction flux(&fes);
      flux.SetFromTrueDofs(f);
      return flux;
    };

    auto smooth_flux_func = build_func(smooth_flux, smooth_flux_fes.GetFinestFESpace());

    auto estimates = ComputeElementLpErrors(smooth_flux_func, normp, coef);

    // Normalize the error by the solution L2 norm, ensures reductions are well scaled.
    const auto normalization = ComputeScalarLpNorm(field, 2);
    std::for_each(estimates.begin(), estimates.end(),
                  [&normalization](auto &x) { x /= normalization; });

    // Mpi::Print("|| u ||_L2 : {}\n", normalization);
    return estimates;
  }
}

namespace
{
mfem::Vector ComputeElementLpErrors(const mfem::ParGridFunction &sol, double p,
                                    mfem::VectorCoefficient &exsol,
                                    int quad_order_increment)
{
  auto &fes = *sol.ParFESpace();

  MFEM_VERIFY(p < mfem::infinity(), "Do not use this routine to compute Linf norm\n");

  mfem::Vector error(fes.GetNE());
  error = 0.0;

  mfem::DenseMatrix vals, exact_vals;
  mfem::Vector loc_errs;

  for (int i = 0; i < fes.GetNE(); ++i)
  {
    const auto &fe = *fes.GetFE(i);
    auto &T = *fes.GetElementTransformation(i);
    const auto &ir = *utils::GetDefaultRule(fe, T, quad_order_increment);

    sol.GetVectorValues(T, ir, vals);
    exsol.Eval(exact_vals, T, ir);

    vals -= exact_vals;

    for (int j = 0; j < ir.GetNPoints(); ++j)
    {
      // Combine components
      double component = 0.0;
      for (int c = 0; c < vals.Height(); ++c)
      {
        component += std::pow(std::abs(vals(c, j)), p);
      }

      const auto &ip = ir.IntPoint(j);
      T.SetIntPoint(&ip);

      error[i] += ip.weight * T.Weight() * component;
    }
    error[i] = std::pow(std::abs(error[i]), 1 / p);
  }
  return error;
}

double ComputeVectorLpNorm(const mfem::ParGridFunction &sol, double p,
                           int quad_order_increment)
{
  auto &fes = *sol.ParFESpace();
  const int nelem = fes.GetNE();

  MFEM_VERIFY(p < mfem::infinity(), "Do not use this routine to compute Linf norm\n");

  mfem::DenseMatrix vals;
  mfem::Vector loc, elem_norm(nelem);

  elem_norm = 0.0;

  for (int i = 0; i < fes.GetNE(); ++i)
  {
    const auto &fe = *fes.GetFE(i);
    auto &T = *fes.GetElementTransformation(i);
    const auto &ir = *utils::GetDefaultRule(fe, T, quad_order_increment);
    sol.GetVectorValues(T, ir, vals);

    for (int j = 0; j < ir.GetNPoints(); ++j)
    {
      // Combine components
      double component = 0.0;
      for (int c = 0; c < vals.Height(); ++c)
      {
        component += std::pow(std::abs(vals(c, j)), p);
      }

      const auto &ip = ir.IntPoint(j);
      T.SetIntPoint(&ip);

      elem_norm[i] += ip.weight * T.Weight() * component;
    }

    // negative quadrature weights might case elemental norms to have been
    // negative, correct this before accumulation.
    elem_norm[i] = std::abs(elem_norm[i]);
  }

  // (sum_e (\| u \|_Lp(e))^p)^(1/p)
  auto norm = std::accumulate(elem_norm.begin(), elem_norm.end(), 0.0);

  Mpi::GlobalSum(1, &norm, Mpi::World());

  norm = std::pow(norm, 1.0 / p);

  return norm;
}

double ComputeScalarLpNorm(const mfem::ParGridFunction &sol, double p,
                           int quad_order_increment)
{
  auto &fes = *sol.ParFESpace();
  const int nelem = fes.GetNE();

  MFEM_VERIFY(p < mfem::infinity(), "Do not use this routine to compute Linf norm\n");

  mfem::Vector vals, elem_norm(nelem);

  elem_norm = 0.0;

  for (int i = 0; i < fes.GetNE(); ++i)
  {
    const auto &fe = *fes.GetFE(i);
    auto &T = *fes.GetElementTransformation(i);
    const auto &ir = *utils::GetDefaultRule(fe, T, quad_order_increment);
    sol.GetValues(T, ir, vals);

    for (int j = 0; j < ir.GetNPoints(); ++j)
    {
      const auto &ip = ir.IntPoint(j);
      T.SetIntPoint(&ip);
      elem_norm[i] += ip.weight * T.Weight() * std::pow(std::abs(vals(j)), p);
    }

    // negative quadrature weights might cause elemental norms to have been
    // negative.
    elem_norm[i] = std::abs(elem_norm[i]);
  }

  // (∑ₑ ∫ₑ|u|^p)^(1/p)
  auto norm = std::accumulate(elem_norm.begin(), elem_norm.end(), 0.0);

  Mpi::GlobalSum(1, &norm, Mpi::World());

  norm = std::pow(norm, 1.0 / p);

  return norm;
}
}  // namespace

}  // namespace palace
