// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "errorindicator.hpp"

#include <mfem/general/forall.hpp>

namespace palace
{

// CUSTOM CONVERGENCE REMOVE
// void ErrorIndicator::AddIndicator(const Vector &indicator)
// {
//   if (n == 0)
//   {
//     local = indicator;
//     n = 1;
//     last_error = std::numeric_limits<double>::max();
//     consecutive_converged = 0;
//     return;
//   }

//   // The average local indicator is used rather than the indicator for the maximum
//   // error to drive the adaptation, to account for a local error that might be marginally
//   // important to many solves, rather than only large in one solve.
//   MFEM_ASSERT(local.Size() == indicator.Size(),
//               "Unexpected size mismatch for ErrorIndicator::AddIndicator!");

//   // The local indicators must be squared before combining, so that the global error
//   // calculation is valid:
//   //                            E = √(1/N ∑ₙ ∑ₖ ηₖₙ²)
//   // from which it follows that:
//   //                            E² = 1/N ∑ₙ ∑ₖ ηₖₙ²
//   //                               = 1/N ∑ₙ Eₙ²
//   // Namely the average of the global error indicators included in the reduction.
//   // Squaring both sides means the summation can be rearranged, and then the local error
//   // indicators become:
//   //                            eₖ = √(1/N ∑ₙ ηₖₙ²)
//   const bool use_dev = local.UseDevice() || indicator.UseDevice();
//   const int N = local.Size();
//   const int Dn = n;
//   const auto *DI = indicator.Read();
//   auto *DL = local.ReadWrite();
//   const double weight = is_jj ? jj_weight : 1.0;
//   // mfem::forall_switch(
//   //     use_dev, N, [=] MFEM_HOST_DEVICE(int i)
//   //     { DL[i] = std::sqrt((DL[i] * DL[i] * Dn + DI[i] * DI[i]) / (Dn + 1)); });
//   // CUSTOM CONVERGENCE
//   mfem::forall_switch(
//   use_dev, N, [=] MFEM_HOST_DEVICE(int i) {
//     const double weighted_ind = DI[i] * weight;
//     DL[i] = std::sqrt((DL[i] * DL[i] * Dn + 
//                       weighted_ind * weighted_ind) / (Dn + 1));
//   });
//   double current_error = linalg::Norml2(GetComm(), local);
//   if (current_error < last_error * (1.0 + rel_tol)) {
//     consecutive_converged++;
//   } else {
//     consecutive_converged = 0;
//   }
//   last_error = current_error;

//   // More samples have been added, update for the running average.
//   n += 1;
// }

// CUSTOM CONVERGENCE
void ErrorIndicator::AddIndicator(const Vector &indicator, bool is_jj) {
  if (n == 0) {
    local = indicator;
    n = 1;
    last_error = std::numeric_limits<double>::max();
    consecutive_converged = 0;
    return;
  }

  // Higher weight for JJ regions
  const double weight = is_jj ? jj_weight : 1.0;

  const bool use_dev = local.UseDevice() || indicator.UseDevice();
  const int N = local.Size();
  const auto *DI = indicator.Read();
  const auto *DL = local.Read();
  
  // Calculate MSE without modifying local yet
  Vector diff_sq(N);
  diff_sq.UseDevice(use_dev);
  auto *DS = diff_sq.Write();
  
  mfem::forall_switch(
    use_dev, N, [=] MFEM_HOST_DEVICE(int i) {
      // Square of difference for MSE with weighting
      const double diff = (DI[i] - DL[i]) * weight;
      DS[i] = diff * diff;
    });

  // Sum up MSE and divide by number of elements
  double current_error = linalg::Sum(MPI_COMM_WORLD, diff_sq) / N;

  // Check for convergence criteria
  if (current_error < global_tol && 
      std::abs(current_error - last_error)/std::abs(last_error) < relative_tol) {
    consecutive_converged++;
  } else {
    consecutive_converged = 0;
  }
  
  // Store the current indicator for next comparison
  local = indicator;
  
  // Update tracking variables
  last_error = current_error;
  n++;
}

}  // namespace palace
