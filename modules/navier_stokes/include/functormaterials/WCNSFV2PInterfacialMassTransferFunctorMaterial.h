//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FunctorMaterial.h"

/**
 * Computes the interfacial mass transfer rate of the two-phase mixture model from the transported
 * interfacial area concentration, rather than taking it as a prescribed rate.
 *
 * The rate follows the interfacial energy jump condition: the heat the interface exchanges with
 * the surrounding fluid is what converts mass from one phase to the other, so
 *
 * \f[
 *   \Gamma = \frac{\chi_p h_i \left(T - T_{sat}\right)}{h_{fg}}
 * \f]
 *
 * in which \f$ \chi_p \f$ is the interfacial area concentration, \f$ h_i \f$ the interfacial heat
 * transfer coefficient per unit area, and \f$ h_{fg} \f$ the latent heat. Positive generates the
 * dispersed phase.
 *
 * The mixture model carries one energy equation and therefore one temperature, so the two phases
 * are in thermal equilibrium with each other and the driving potential is the departure of that
 * shared temperature from saturation rather than a difference between phase temperatures. The
 * interface itself is at saturation, so the enthalpy the transfer has to supply is \f$ h_{fg} \f$.
 * This is the same driving potential the mixture-model phase change closures of the CFD codes use,
 * with the rate coefficient computed from the solved area instead of being tuned.
 *
 * The coefficient is the bubbly-flow interfacial Nusselt number of RELAP5/MOD3, the modified
 * Lee-Ryley correlation of Section 4.1.1.1.1 of NUREG/CR-5535 Volume 4,
 *
 * \f[
 *   h_i = \frac{k_c}{d_b}\left(2 + 0.74\, Re_b^{1/2}\right) ,
 *   \qquad d_b = \frac{\psi \alpha}{\chi_p} ,
 *   \qquad Re_b = \frac{\rho_c d_b \left|u_s\right|}{\mu_c} ,
 * \f]
 *
 * with the Prandtl number dependence dropped, as that reference drops it for bubbly flow. RELAP5
 * forms its volumetric coefficient as the product of this with an interfacial area of
 * \f$ 3.6\alpha/d_b \f$ obtained from a critical Weber number. Here the transported area takes the
 * place of that algebraic estimate, which is the one substitution this closure makes: the heat
 * transfer physics is the reference correlation, the geometry is solved for.
 *
 * Two limitations follow from that choice and are worth stating. The correlation was assessed
 * against RELAP5's own area, so the product is no longer the quantity that reference validated.
 * And RELAP5 takes the larger of this and a Plesset-Zwick bubble growth rate, which dominates at
 * strong superheat; only the Lee-Ryley branch is evaluated here, so rapid flashing is outside what
 * this closure has been written for.
 *
 * Wall nucleation is not included. This is bulk transfer between phases already in contact, which
 * grows the particles present at fixed number density and is what the two thirds exponent of the
 * interfacial area source assumes. Vapour generated at a heated wall creates new particles instead
 * and needs a nucleation source the area equation does not have.
 */
class WCNSFV2PInterfacialMassTransferFunctorMaterial : public FunctorMaterial
{
public:
  static InputParameters validParams();

  WCNSFV2PInterfacialMassTransferFunctorMaterial(const InputParameters & parameters);

protected:
  /// The dimension of the simulation
  const unsigned int _dim;

  /// Shape factor of the averaged particle size, 6 for spheres
  const Real _shape_factor;

  /// Interfacial area concentration, the transported area this closure is built on
  const Moose::Functor<Real> & _interfacial_area;

  /// Volume fraction of the dispersed phase
  const Moose::Functor<Real> & _fd;

  /// Temperature of the mixture, shared by the phases
  const Moose::Functor<Real> & _temperature;

  /// Saturation temperature of the transition the transfer represents
  const Moose::Functor<Real> & _T_saturation;

  /// Latent heat of that transition, per unit mass
  const Moose::Functor<Real> & _latent_heat;

  /// Continuous phase properties, which carry the interfacial heat transfer
  const Moose::Functor<Real> & _rho_c;
  const Moose::Functor<Real> & _mu_c;
  const Moose::Functor<Real> & _k_c;

  /// Components of the slip velocity, which set the particle Reynolds number
  const Moose::Functor<Real> & _u_slip;
  const Moose::Functor<Real> * const _v_slip;
  const Moose::Functor<Real> * const _w_slip;
};
