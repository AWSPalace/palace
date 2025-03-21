// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_ERROR_INDICATORS_HPP
#define PALACE_FEM_ERROR_INDICATORS_HPP

#include <array>
#include <vector>
#include "linalg/vector.hpp"
#include "utils/communication.hpp"

namespace palace
{

//
// Storage for error estimation results from a simulation which involves one or more solves,
// required in the AMR loop.
//
class ErrorIndicator
{
protected:
  // Elemental localized error indicators. Used for marking elements for
  // refinement and coarsening.
  Vector local;

  // Number of samples.
  int n;

  // CUSTOM CONVERGENCE
  int consecutive_converged{0};
  int required_consecutive{3}; 
  double last_error{std::numeric_limits<double>::max()};
  double global_tol{1e-2};
  double relative_tol{1e-3};
  double jj_weight{10.0};  // Significantly increased weight for JJ elements

public:
  ErrorIndicator(Vector &&local) : local(std::move(local)), n(1)
  {
    this->local.UseDevice(true);
  }
  ErrorIndicator() : n(0) { local.UseDevice(true); }

  // Add an indicator to the running total.
  // CUSTOM CONVERGENCE REMOVE
  //void AddIndicator(const Vector &indicator);

  // Return the local error indicator.
  const auto &Local() const { return local; }

  // Return the global error indicator.
  auto Norml2(MPI_Comm comm) const { return linalg::Norml2(comm, local); }

  // Return the largest local error indicator.
  auto Max(MPI_Comm comm) const
  {
    auto max = local.Max();
    Mpi::GlobalMax(1, &max, comm);
    return max;
  }

  // Return the smallest local error indicator.
  auto Min(MPI_Comm comm) const
  {
    auto min = local.Min();
    Mpi::GlobalMin(1, &min, comm);
    return min;
  }

  // Return the mean local error indicator.
  auto Mean(MPI_Comm comm) const { return linalg::Mean(comm, local); }
  // CUSTOM CONVERGENCE
  void SetConvergenceParams(double tol, double rel_tol, int required_consec) {
    global_tol = tol;
    relative_tol = rel_tol;
    required_consecutive = required_consec;
  }

  bool HasConverged() const {
    return consecutive_converged >= required_consecutive;
  }

  // Modify signature of existing AddIndicator:
  void AddIndicator(const Vector &indicator, bool is_jj = false);
};

}  // namespace palace

#endif  // PALACE_FEM_ERROR_INDICATORS_HPP
