// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "lumpedportoperator.hpp"

#include "fem/coefficient.hpp"
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

namespace palace
{

using namespace std::complex_literals;

LumpedPortData::LumpedPortData(const config::LumpedPortData &data,
                               const MaterialOperator &mat_op, const mfem::ParMesh &mesh)
  : mat_op(mat_op)
{
  // Check inputs. Only one of the circuit or per square properties should be specified
  // for the port boundary.
  bool has_circ = (std::abs(data.R) + std::abs(data.L) + std::abs(data.C) > 0.0);
  bool has_surf = (std::abs(data.Rs) + std::abs(data.Ls) + std::abs(data.Cs) > 0.0);
  MFEM_VERIFY(has_circ || has_surf,
              "Lumped port boundary has no R/L/C or Rs/Ls/Cs defined, needs "
              "at least one!");
  MFEM_VERIFY(!(has_circ && has_surf),
              "Lumped port boundary has both R/L/C and Rs/Ls/Cs defined, "
              "should only use one!");
  excitation = data.excitation;
  active = data.active;
  if (excitation)
  {
    if (has_circ)
    {
      MFEM_VERIFY(data.R > 0.0, "Excited lumped port must have nonzero resistance!");
      MFEM_VERIFY(data.C == 0.0 && data.L == 0.0,
                  "Lumped port excitations do not support nonzero reactance!");
    }
    else
    {
      MFEM_VERIFY(data.Rs > 0.0, "Excited lumped port must have nonzero resistance!");
      MFEM_VERIFY(data.Cs == 0.0 && data.Ls == 0.0,
                  "Lumped port excitations do not support nonzero reactance!");
    }
  }

  // Construct the lumped element geometry information for the port, converting circuit
  // parameters like R to per-square values Rs if circuit values were specified.
  if (data.mode.has_value())
  {
    // Coaxial mode (radial direction).
    elems.push_back(std::make_unique<CoaxialElementData>(data.mode.value(), data.attributes,
                                                         mesh));
    if (has_circ)
    {
      // Conversion factor: Z_sq = R * (2π / ln(r_o / r_i)).
      const double f = 2.0 * M_PI / elems.back()->GetGeometryLength();  // [1]
      R = data.R * f;                                                    // [Ω/sq]
      L = data.L * f;                                                    // [H/sq]
      C = data.C / f;                                                    // [F·sq]
    }
    else
    {
      R = data.Rs;  // [Ω/sq]
      L = data.Ls;  // [H/sq]
      C = data.Cs;  // [F·sq]
    }
  }
  else
  {
    // Direction-based mode.
    std::array<double, 3> dir = {data.direction[0], data.direction[1], data.direction[2]};
    elems.push_back(
        std::make_unique<UniformElementData>(dir, data.attributes, mesh));
    if (has_circ)
    {
      // Conversion factor: Z_sq = R * (w / l).
      const double f = GetToSquare(*elems.back());  // [1]
      R = data.R * f;                               // [Ω/sq]
      L = data.L * f;                               // [H/sq]
      C = data.C / f;                               // [F·sq]
    }
    else
    {
      R = data.Rs;  // [Ω/sq]
      L = data.Ls;  // [H/sq]
      C = data.Cs;  // [F·sq]
    }
  }

  // For the case of a multi-element port, add additional port elements if they are
  // specified.
  if (!data.elements.empty())
  {
    for (const auto &elem : data.elements)
    {
      std::array<double, 3> dir = {elem.direction[0], elem.direction[1], elem.direction[2]};
      elems.push_back(
          std::make_unique<UniformElementData>(dir, elem.attributes, mesh));
    }
  }
}

void LumpedPortData::InitializeLinearForms(mfem::ParFiniteElementSpace &nd_fespace) const
{
  if (s && v)
  {
    return;
  }
  s = std::make_unique<mfem::LinearForm>(&nd_fespace);
  v = std::make_unique<mfem::LinearForm>(&nd_fespace);
  for (const auto &elem : elems)
  {
    for (int attr : elem->GetAttrList())
    {
      auto s_coeff = elem->GetModeCoefficient(1.0);
      s->AddBdrFaceIntegrator(new SurfaceVectorIntegrator(*s_coeff), attr);
      auto v_coeff = elem->GetModeCoefficient(1.0 / elems.size());
      v->AddBdrFaceIntegrator(new SurfaceVectorIntegrator(*v_coeff), attr);
    }
  }
  s->Assemble();
  v->Assemble();
}

std::complex<double>
LumpedPortData::GetCharacteristicImpedance(double omega, Branch branch) const
{
  std::complex<double> Z = 0.0;
  switch (branch)
  {
    case Branch::TOTAL:
      if (std::abs(omega) < 1.0e-15)
      {
        // DC case.
        Z = R;
      }
      else
      {
        if (std::abs(L) < 1.0e-15 && std::abs(C) < 1.0e-15)
        {
          // Resistor.
          Z = R;
        }
        else if (std::abs(R) < 1.0e-15 && std::abs(C) < 1.0e-15)
        {
          // Inductor.
          Z = 1i * omega * L;
        }
        else if (std::abs(R) < 1.0e-15 && std::abs(L) < 1.0e-15)
        {
          // Capacitor.
          Z = 1.0 / (1i * omega * C);
        }
        else if (std::abs(C) < 1.0e-15)
        {
          // Resistor + inductor in series.
          Z = R + 1i * omega * L;
        }
        else if (std::abs(L) < 1.0e-15)
        {
          // Resistor + capacitor in series.
          Z = R + 1.0 / (1i * omega * C);
        }
        else if (std::abs(R) < 1.0e-15)
        {
          // Inductor + capacitor in series (resonator).
          Z = 1i * omega * L + 1.0 / (1i * omega * C);
        }
        else
        {
          // Resistor + inductor + capacitor in series (lossy resonator).
          Z = R + 1i * omega * L + 1.0 / (1i * omega * C);
        }
      }
      break;
    case Branch::R:
      Z = R;
      break;
    case Branch::L:
      Z = (std::abs(omega) > 1.0e-15) ? 1i * omega * L : 0.0;
      break;
    case Branch::C:
      Z = (std::abs(omega) > 1.0e-15) ? 1.0 / (1i * omega * C) : 0.0;
      break;
  }
  return Z;
}

double LumpedPortData::GetExcitationPower() const
{
  return 0.5;  // Power normalization constant for excitation.
}

double LumpedPortData::GetExcitationVoltage() const
{
  return 1.0;  // Voltage normalization constant for excitation.
}

std::complex<double> LumpedPortData::GetPower(GridFunction &E, GridFunction &B) const
{
  // Compute the complex power flow through the port.
  InitializeLinearForms(*E.ParFESpace());
  std::complex<double> dot(0.0, 0.0);
  if (std::abs(R) > 0.0)
  {
    // S = 0.5 * V^2 / R.
    const std::complex<double> V_i = GetVoltage(E);
    dot = 0.5 * V_i * std::conj(V_i) / R;
  }
  else if (std::abs(L) > 0.0)
  {
    // S = 0.5 * omega * I^2 * L.
    const std::complex<double> I_i = (*s) * E.Real();
    dot = 0.5 * I_i * std::conj(I_i) * L;
  }
  else if (std::abs(C) > 0.0)
  {
    // S = 0.5 * V^2 / (omega * C).
    const std::complex<double> V_i = GetVoltage(E);
    dot = 0.5 * V_i * std::conj(V_i) / C;
  }
  Mpi::GlobalSum(1, &dot, E.GetComm());
  return dot;
}

std::complex<double> LumpedPortData::GetSParameter(GridFunction &E) const
{
  // Compute the reflection parameter at the port.
  InitializeLinearForms(*E.ParFESpace());
  std::complex<double> dot((*s) * E.Real(), 0.0);
  if (E.HasImag())
  {
    dot.imag((*s) * E.Imag());
  }
  Mpi::GlobalSum(1, &dot, E.GetComm());
  return dot;
}

std::complex<double> LumpedPortData::GetVoltage(GridFunction &E) const
{
  // Compute the average voltage across the port.
  InitializeLinearForms(*E.ParFESpace());
  std::complex<double> dot((*v) * E.Real(), 0.0);
  if (E.HasImag())
  {
    dot.imag((*v) * E.Imag());
  }
  Mpi::GlobalSum(1, &dot, E.GetComm());
  return dot;
}

double LumpedPortData::ComputeElectricFieldEnergy(const mfem::Vector &field_mag) const
{
  // Since we don't have direct access to element indices, we'll use a simple
  // approximation for junction energy by checking if the port contains attributes
  // that match those in field_mag
  
  // In a more advanced implementation, we would integrate the field energy
  // over the junction volume, but here we'll just approximate
  
  double total_energy = 0.0;
  
  // Simplistic approach - if this is an inductive element (L > 0)
  // add energy contribution proportional to L
  if (std::abs(L) > 0.0) {
    total_energy = L * 1e-3; // Arbitrary scaling factor 
  }
  
  return total_energy;
}

LumpedPortOperator::LumpedPortOperator(const IoData &iodata, const MaterialOperator &mat_op,
                                       const mfem::ParMesh &mesh)
{
  // Set up lumped port boundary conditions.
  SetUpBoundaryProperties(iodata, mat_op, mesh);
  PrintBoundaryInfo(iodata, mesh);
}

void LumpedPortOperator::SetUpBoundaryProperties(const IoData &iodata,
                                                 const MaterialOperator &mat_op,
                                                 const mfem::ParMesh &mesh)
{
  // Check that lumped port boundary attributes have been specified correctly.
  if (!iodata.boundaries.lumpedport.empty())
  {
    int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
    mfem::Array<int> bdr_attr_marker(bdr_attr_max), port_marker(bdr_attr_max);
    bdr_attr_marker = 0;
    port_marker = 0;
    for (auto attr : mesh.bdr_attributes)
    {
      MFEM_ASSERT(0 < attr && attr <= bdr_attr_max, "Invalid boundary attribute in mesh!");
      bdr_attr_marker[attr - 1] = 1;
    }
    int offset = 0;
    for (const auto &data : iodata.boundaries.lumpedport)
    {
      // Check for reasonable index number to identify the port.
      MFEM_VERIFY(data.idx > 0, "Lumped port has non-positive index {:d}!", data.idx);
      MFEM_VERIFY(data.attributes.size() == 1 || data.elements.size() == 0,
                  "For port [{:d}], to use multiple port elements, remove the "
                  "\"Attributes\" field and use only \"Elements\"!",
                  data.idx);

      // Check for duplicate port index numbers.
      MFEM_VERIFY(port_op.count(data.idx) == 0, "Duplicate lumped port index {:d}!", data.idx);

      // For non-empty port elements, check reasonable attribute list.
      if (!data.elements.empty())
      {
        for (const auto &elem : data.elements)
        {
          const int elen = static_cast<int>(elem.attributes.size());
          MFEM_VERIFY(elen > 0,
                      "Invalid empty attribute list for a port [{:d}] element!", data.idx);
          // Check for missing or invalid attributes.
          for (auto a : elem.attributes)
          {
            MFEM_VERIFY(a > 0 && a <= bdr_attr_max,
                        "Invalid lumped port attribute {:d} for port [{:d}]!", a, data.idx);
            MFEM_VERIFY(bdr_attr_marker[a - 1] == 1,
                        "Missing lumped port attribute {:d} in mesh for port [{:d}]!",
                        a, data.idx);
          }
          for (size_t i = 0; i < elem.attributes.size(); i++)
          {
            // Check for reused attributes.
            for (const auto &[idx, p] : port_op)
            {
              const auto &pattr = p.elems.front()->GetAttrList();
              for (int k = 0; k < pattr.Size(); k++)
              {
                for (auto a : elem.attributes)
                {
                  MFEM_VERIFY(pattr[k] != a,
                              "Attribute {:d} reused between port [{:d}] and port [{:d}]!",
                              a, p.idx, data.idx);
                }
              }
            }
            port_marker[elem.attributes[i] - 1] = 1;
          }
        }
      }
      else
      {
        const int alen = static_cast<int>(data.attributes.size());
        MFEM_VERIFY(alen > 0,
                    "Invalid empty attribute list for a lumped port boundary!");
        // Check for missing or invalid attributes.
        for (auto a : data.attributes)
        {
          MFEM_VERIFY(a > 0 && a <= bdr_attr_max,
                      "Invalid lumped port attribute {:d} for port [{:d}]!", a, data.idx);
          MFEM_VERIFY(bdr_attr_marker[a - 1] == 1,
                      "Missing lumped port attribute {:d} in mesh for port [{:d}]!",
                      a, data.idx);
        }
        // Check for reused attributes.
        for (const auto &[idx, p] : port_op)
        {
          const auto &pattr = p.elems.front()->GetAttrList();
          for (int k = 0; k < pattr.Size(); k++)
          {
            for (auto a : data.attributes)
            {
              MFEM_VERIFY(pattr[k] != a,
                          "Attribute {:d} reused between port [{:d}] and port [{:d}]!",
                          a, idx, data.idx);
            }
          }
        }

        for (auto a : data.attributes)
        {
          port_marker[a - 1] = 1;
        }
      }

      // Add the port operators to the list.
      const auto &it = port_op.insert({data.idx, {data, mat_op, mesh}});
      if (data.excitation)
      {
        // The first half of the two-element excitation source array corresponds to the
        // port's associated lumped element. The second half corresponds to the
        // actual solution variable we're solving for in the system.
        source.resize(std::max(static_cast<int>(source.size()), 2 * (offset + 1)));
        source[2 * offset] = &it.first->second;
        source[2 * offset + 1] = &it.first->second;
        offset += 1;
      }
    }
  }
}

void LumpedPortOperator::PrintBoundaryInfo(const IoData &iodata,
                                           const mfem::ParMesh &mesh) const
{
  if (!port_op.empty())
  {
    if (mesh.GetMyRank() == 0)
    {
      if (port_op.size() == 1)
      {
        const auto &data = port_op.begin()->second;
        if (data.excitation)
        {
          Mpi::Print("\nExcited lumped port boundary:\n");
        }
        else
        {
          Mpi::Print("\nLumped port boundary:\n");
        }
      }
      else
      {
        bool has_exc = false;
        for (const auto &port : port_op)
        {
          has_exc = has_exc || port.second.excitation;
        }
        if (has_exc)
        {
          Mpi::Print("\nExcited lumped port boundaries ({:d}):\n", port_op.size());
        }
        else
        {
          Mpi::Print("\nLumped port boundaries ({:d}):\n", port_op.size());
        }
      }
      for (const auto &[idx, data] : port_op)
      {
        bool is_active = data.active;
        bool is_excited = data.excitation;
        bool is_res = std::abs(data.R) > 0.0 && std::abs(data.L) > 0.0 && std::abs(data.C) > 0.0;
        bool is_cap = std::abs(data.R) <= 0.0 && std::abs(data.L) <= 0.0 && std::abs(data.C) > 0.0;
        bool is_ind = std::abs(data.R) <= 0.0 && std::abs(data.L) > 0.0 && std::abs(data.C) <= 0.0;
        bool is_res_lc = std::abs(data.R) <= 0.0 && std::abs(data.L) > 0.0 && std::abs(data.C) > 0.0;
        bool is_res_lr = std::abs(data.R) > 0.0 && std::abs(data.L) > 0.0 && std::abs(data.C) <= 0.0;
        bool is_res_rc = std::abs(data.R) > 0.0 && std::abs(data.L) <= 0.0 && std::abs(data.C) > 0.0;
        bool is_res_r = std::abs(data.R) > 0.0 && std::abs(data.L) <= 0.0 && std::abs(data.C) <= 0.0;
        const double f_GHz = 1.0e-9 * iodata.GetFrequencyPointOmega() / (2.0 * M_PI);
        const double abs_omega = std::abs(iodata.GetFrequencyPointOmega());

        std::string desc;
        if (is_active)
        {
          desc = is_excited ? "Excited" : "Passive";
        }
        else
        {
          desc = "Deactivated";
        }
        if (is_res)
        {
          const double f_res = 1.0e-9 / (2.0 * M_PI * std::sqrt(data.L * data.C));
          Mpi::Print(" [{:d}]  {:<11s}  resistor + inductor + capacitor (lossy resonator):\n"
                     "                 R = {:.4g} Ω/sq,  L = {:.4g} H/sq,  C = {:.4g} F·sq\n"
                     "                 f_res = {:.4g} GHz,  Q = {:.4g}\n",
                     idx, desc,
                     data.R, data.L, data.C,
                     f_res, 1.0 / data.R * std::sqrt(data.L / data.C));
          if (abs_omega > 0.0)
          {
            const double omega_res = 1.0 / std::sqrt(data.L * data.C);
            const double delta_f = f_GHz - (1.0e-9 * omega_res / (2.0 * M_PI));
            double Z = std::abs(data.GetCharacteristicImpedance(iodata.GetFrequencyPointOmega()));
            Mpi::Print("                 (f = {:.4g} GHz, Δf/f_res = {:.2g},  |Z| = {:.4g} Ω/sq)\n",
                       f_GHz, delta_f / f_res, Z);
          }
        }
        else if (is_cap)
        {
          Mpi::Print(" [{:d}]  {:<11s}  capacitor:\n"
                     "                 C = {:.4g} F·sq\n",
                     idx, desc,
                     data.C);
          if (abs_omega > 0.0)
          {
            const double Z = 1.0 / (abs_omega * data.C);
            Mpi::Print("                 (f = {:.4g} GHz,  |Z| = {:.4g} Ω/sq)\n",
                       f_GHz, Z);
          }
        }
        else if (is_ind)
        {
          Mpi::Print(" [{:d}]  {:<11s}  inductor:\n"
                     "                 L = {:.4g} H/sq\n",
                     idx, desc,
                     data.L);
          if (abs_omega > 0.0)
          {
            const double Z = abs_omega * data.L;
            Mpi::Print("                 (f = {:.4g} GHz,  |Z| = {:.4g} Ω/sq)\n",
                       f_GHz, Z);
          }
        }
        else if (is_res_lc)
        {
          const double f_res = 1.0e-9 / (2.0 * M_PI * std::sqrt(data.L * data.C));
          Mpi::Print(" [{:d}]  {:<11s}  inductor + capacitor (resonator):\n"
                     "                 L = {:.4g} H/sq,  C = {:.4g} F·sq\n"
                     "                 f_res = {:.4g} GHz\n",
                     idx, desc,
                     data.L, data.C,
                     f_res);
          if (abs_omega > 0.0)
          {
            const double omega_res = 1.0 / std::sqrt(data.L * data.C);
            const double delta_f = f_GHz - (1.0e-9 * omega_res / (2.0 * M_PI));
            const double Z = std::abs(data.GetCharacteristicImpedance(iodata.GetFrequencyPointOmega()));
            Mpi::Print("                 (f = {:.4g} GHz, Δf/f_res = {:.2g},  |Z| = {:.4g} Ω/sq)\n",
                       f_GHz, delta_f / f_res, Z);
          }
        }
        else if (is_res_lr)
        {
          Mpi::Print(" [{:d}]  {:<11s}  resistor + inductor:\n"
                     "                 R = {:.4g} Ω/sq,  L = {:.4g} H/sq\n",
                     idx, desc,
                     data.R, data.L);
          if (abs_omega > 0.0)
          {
            const double Z = std::abs(data.GetCharacteristicImpedance(iodata.GetFrequencyPointOmega()));
            Mpi::Print("                 (f = {:.4g} GHz,  |Z| = {:.4g} Ω/sq)\n",
                       f_GHz, Z);
          }
        }
        else if (is_res_rc)
        {
          Mpi::Print(" [{:d}]  {:<11s}  resistor + capacitor:\n"
                     "                 R = {:.4g} Ω/sq,  C = {:.4g} F·sq\n",
                     idx, desc,
                     data.R, data.C);
          if (abs_omega > 0.0)
          {
            const double Z = std::abs(data.GetCharacteristicImpedance(iodata.GetFrequencyPointOmega()));
            Mpi::Print("                 (f = {:.4g} GHz,  |Z| = {:.4g} Ω/sq)\n",
                       f_GHz, Z);
          }
        }
        else if (is_res_r)
        {
          Mpi::Print(" [{:d}]  {:<11s}  resistor:\n"
                     "                 R = {:.4g} Ω/sq\n",
                     idx, desc,
                     data.R);
        }
      }
    }
  }
}

bool LumpedPortOperator::HasLumpedPorts() const
{
  return !port_op.empty();
}

void LumpedPortOperator::ApplyEssentialBoundaryConditions(
    mfem::SparseMatrix &mat, mfem::Vector &rhs,
    const mfem::Array<int> &tdof_list) const
{
  // Apply the boundary conditions for the lumped resistors.
  for (const auto &[idx, data] : port_op)
  {
    if (std::abs(data.R) > 0.0 && data.active)
    {
      for (const auto &elem : data.elems)
      {
        const auto &attr_list = elem->GetAttrList();
        for (int i = 0; i < attr_list.Size(); i++)
        {
          mat.EliminateBCFromDofs(tdof_list, attr_list[i], rhs);
        }
      }
    }
  }
}

}  // namespace palace