//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "LinearFVElementalKernel.h"

/**
 * Sources and sinks of the one-group interfacial area concentration transport equation, for the
 * linear finite volume discretization.
 *
 * The transport equation, in conservative form,
 *
 * \f[
 *   \frac{\partial (\rho_g \chi_p)}{\partial t}
 *     + \nabla \cdot \left( \rho_g \mathbf{u}_g \chi_p \right)
 *   = \frac{1}{3}\frac{D\rho_g}{Dt}\chi_p
 *     + \frac{2}{3}\frac{\dot m_g}{\alpha_g}\chi_p
 *     + \rho_g \left( S_{RC} + S_{WE} + S_{TI} \right) ,
 * \f]
 *
 * in which \f$ \chi_p \f$ is the interfacial area concentration, \f$ \rho_g \f$ and
 * \f$ \alpha_g \f$ the density and volume fraction of the dispersed (gas) phase, and
 * \f$ \dot m_g \f$ the mass transfer rate into the gas phase per unit mixture volume. The first two
 * terms on the right are the change of area from expansion of the dispersed phase and from phase
 * change. \f$ S_{RC} \f$, \f$ S_{WE} \f$ and \f$ S_{TI} \f$ are the coalescence sinks from random
 * collision and wake entrainment, and the breakage source from turbulent impact.
 *
 * This object assembles everything except the time derivative and the advection, which are the
 * business of LinearFVTimeDerivative given a factor of \f$ \rho_g \f$ and of
 * LinearFVScalarAdvection given a density of \f$ \rho_g \f$ and the dispersed phase velocity.
 *
 * Two closure sets are offered, selected by the 'model' parameter. The averaged particle size is
 * common to both, \f$ d_b = \psi \alpha_g / \chi_p \f$.
 *
 * **Hibiki and Ishii**, Int. J. Heat Mass Transfer 43 (2000) 2711. There is no wake
 * entrainment model, so
 * \f$ S_{WE} = 0 \f$.
 *
 * \f[
 *   S_{RC} = -\left(\frac{\alpha_g}{\chi_p}\right)^2
 *     \frac{\Gamma_C \alpha_g^2 \epsilon^{1/3}}{d_b^{11/3}\left(\alpha_{g,max}-\alpha_g\right)}
 *     \exp\left(-K_C \frac{d_b^{5/6}\rho_f^{1/2}\epsilon^{1/3}}{\sigma^{1/2}}\right)
 * \f]
 * \f[
 *   S_{TI} = \left(\frac{\alpha_g}{\chi_p}\right)^2
 *     \frac{\Gamma_B \alpha_g\left(1-\alpha_g\right)\epsilon^{1/3}}
 *          {d_b^{11/3}\left(\alpha_{g,max}-\alpha_g\right)}
 *     \exp\left(-K_B \frac{\sigma}{\rho_f d_b^{5/3}\epsilon^{2/3}}\right)
 * \f]
 *
 * **Ishii and Kim**, with the mean bubble fluctuating velocity
 * \f$ u_t = \epsilon^{1/3} d_b^{1/3} \f$, the terminal velocity \f$ u_r \f$ and the Weber number
 * \f$ We = \rho_f u_t^2 d_b / \sigma \f$. The breakage rate is zero below the critical Weber
 * number.
 *
 * \f[
 *   S_{RC} = -\frac{1}{3\pi} C_{RC} u_t \chi_p^2
 *     \left[\frac{1}{\alpha_{g,max}^{1/3}\left(\alpha_{g,max}^{1/3}-\alpha_g^{1/3}\right)}\right]
 *     \left[1 - \exp\left(-C\frac{\alpha_{g,max}^{1/3}\alpha_g^{1/3}}
 *                                {\alpha_{g,max}^{1/3}-\alpha_g^{1/3}}\right)\right]
 * \f]
 * \f[
 *   S_{WE} = -\frac{1}{3\pi} C_{WE} u_r \chi_p^2 C_D^{1/3} ,
 *   \qquad
 *   S_{TI} = \frac{1}{18} C_{TI} u_t \frac{\chi_p^2}{\alpha_g}
 *            \left(1-\frac{We_{cr}}{We}\right)^{1/2}\exp\left(-\frac{We_{cr}}{We}\right)
 * \f]
 */
class LinearWCNSFV2PInterfaceAreaSourceSink : public LinearFVElementalKernel
{
public:
  static InputParameters validParams();
  LinearWCNSFV2PInterfaceAreaSourceSink(const InputParameters & params);

  virtual Real computeMatrixContribution() override;
  virtual Real computeRightHandSideContribution() override;
  virtual void setCurrentElemInfo(const ElemInfo * elem_info) override;

  /**
   * The terminal velocity of a particle, solved together with its drag coefficient.
   *
   * \f$ u_r = \left(d_b g \Delta\rho / (3 C_D \rho_f)\right)^{1/2} \f$ with
   * \f$ C_D = 24(1 + 0.1 Re_D^{0.75})/Re_D \f$ and \f$ Re_D = \rho_f u_r d_b (1-\alpha_g)/\mu_f
   * \f$, so the drag coefficient depends on the velocity it determines. Eliminating \f$ C_D \f$
   * leaves the scalar equation
   *
   * \f[
   *   u_r \left(1 + 0.1 \left(B u_r\right)^{3/4}\right) = T ,
   *   \qquad B = \frac{\rho_f d_b (1-\alpha_g)}{\mu_f} ,
   *   \qquad T = \frac{B}{24}\,\frac{d_b g \Delta\rho}{3 \rho_f} ,
   * \f]
   *
   * whose left hand side is zero at the origin and strictly increasing, so the root is unique and,
   * the bracket factor being at least one, lies in \f$ [0, T] \f$. Solved by a Newton iteration
   * safeguarded by that bracket, in the same way as the algebraic slip closure.
   *
   * @param stokes_velocity the velocity \f$ T \f$ obtained by ignoring the Reynolds correction
   * @param reynolds_per_velocity the factor \f$ B \f$ converting a velocity into \f$ Re_D \f$
   */
  static Real solveTerminalVelocity(Real stokes_velocity, Real reynolds_per_velocity);

protected:
  /// The closure sets offered by the reference
  enum class ModelEnum
  {
    HIBIKI_ISHII = 0,
    ISHII_KIM = 1
  };

  /// Which closure set to evaluate
  const ModelEnum _model;

  /// The dimension of the simulation
  const unsigned int _dim;

  /// Velocity of the dispersed phase, used for the material derivative of its density
  const Moose::Functor<Real> & _u_var;
  const Moose::Functor<Real> * const _v_var;
  const Moose::Functor<Real> * const _w_var;

  /// Density of the dispersed phase, rho_g
  const Moose::Functor<Real> & _rho_d;
  /// Density of the continuous phase, rho_f
  const Moose::Functor<Real> & _rho_f;
  /// Dynamic viscosity of the continuous phase, only needed by the Ishii and Kim drag coefficient
  const Moose::Functor<Real> * const _mu_f;
  /// Volume fraction of the dispersed phase, alpha_g
  const Moose::Functor<Real> & _f_d;
  /// Surface tension between the phases
  const Moose::Functor<Real> & _sigma;
  /// Turbulent dissipation rate of the continuous phase
  const Moose::Functor<Real> & _epsilon;
  /// Mass transfer rate into the dispersed phase per unit mixture volume
  const Moose::Functor<Real> & _mass_transfer_rate;

  /// Shape factor relating the particle size to the phase fraction and the area, psi
  const Real _shape_factor;
  /// Maximum volume fraction admitted by the model
  const Real _f_d_max;
  /// Gravity vector, only needed by the Ishii and Kim terminal velocity
  const RealVectorValue _gravity;

  /// Hibiki and Ishii closure coefficients
  const Real _gamma_c, _kc, _gamma_b, _kb;
  /// Ishii and Kim closure coefficients
  const Real _c_rc, _c_we, _c_ti, _c, _we_cr;

private:
  /// Evaluates both closure sets, filling the coefficients below
  void computeCoefficients();

  /// Coefficient of the terms that multiply the interfacial area concentration directly
  Real _implicit_coefficient;

  /// The remainder, evaluated from the previous iterate
  Real _lagged_source;
};
