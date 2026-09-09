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
#include "MooseLinearVariableFV.h"
#include "NonADFunctorInterface.h"

#include <algorithm>

class Function;

/**
 * Computes the slip velocity of the dispersed phase relative to the continuous phase for the
 * two-phase mixture model, for the linear finite volume discretization, and optionally the
 * diffusion (drift) velocity derived from it.
 *
 * Two distinct relative velocities appear in the mixture model and are easy to confuse. The slip,
 * or relative, velocity is the velocity of the dispersed phase with respect to the continuous
 * phase, \f$ u_s = u_d - u_c \f$, and it is what the algebraic closure of Manninen et al.
 * predicts. The diffusion, or drift, velocity is the velocity of the dispersed phase with respect
 * to the centre of mass of the mixture, \f$ u_{Md} = u_d - u_m \f$, and it is what enters the
 * conservation equations. The two are related by the dispersed phase mass fraction
 * \f$ c_d = \alpha \rho_d / \rho_m \f$,
 *
 * \f[
 *   u_{Md} = \left( 1 - c_d \right) u_s
 * \f]
 *
 * see Manninen, Taivassalo and Kallio, VTT Publications 288 (1996), equation (28). Advecting the
 * phase fraction with \f$ u_m + u_s \f$ rather than \f$ u_m + u_{Md} \f$ is the dilute
 * approximation \f$ c_d \to 0 \f$.
 *
 * This is the linear finite volume counterpart of WCNSFV2PSlipVelocityFunctorMaterial. The two are
 * kept separate rather than sharing one object which switches on the velocity variable type at
 * run time: holding MooseLinearVariableFVReal directly makes the requirements of this closure,
 * cell gradients and a time derivative, a property of the type rather than something rediscovered
 * by a cast on every use.
 */
class LinearWCNSFV2PSlipVelocityFunctorMaterial : public FunctorMaterial,
                                                  public NonADFunctorInterface
{
public:
  static InputParameters validParams();
  LinearWCNSFV2PSlipVelocityFunctorMaterial(const InputParameters & parameters);

protected:
  /// Which drag law closes the force balance
  enum class DragModelEnum
  {
    /// Rigid sphere: Schiller and Naumann below the transition, Newton above
    SCHILLER_NAUMANN = 0,
    /// Distorted fluid particle: C_D grows with size, terminal velocity independent of it
    DISTORTED_PARTICLE = 1,
    /// Whichever of the two resists more, which selects the regime automatically
    AUTOMATIC = 2,
    /// The multi-bubble distorted particle relation of Ishii and Zuber (1979), as written in
    /// equation (57) of Hibiki and Ishii (2003); both papers are cited in full on
    /// ishiiZuberSlipSpeed below
    ISHII_ZUBER = 3
  };

  /// Retrieve a velocity variable and check that it is a linear finite volume variable
  MooseLinearVariableFVReal & getVelocityVariable(const std::string & param_name);

  /**
   * The magnitude of the slip velocity, obtained from the force balance with the drag evaluated
   * at the correct particle Reynolds number.
   *
   * The closure is \f$ u_s = \tau_d / f(Re_p) \cdot K \cdot a \f$ with
   * \f$ Re_p = \rho_c d_d |u_s| / \mu_c \f$, so the drag depends on the very quantity it
   * determines. Writing \f$ s = |u_s| \f$, \f$ R = \rho_c d_d / \mu_c \f$ and
   * \f$ s_0 = \tau_d |K| |a| \f$ for the slip the Stokes limit \f$ f \equiv 1 \f$ would give, the
   * balance is the scalar equation
   *
   * \f[
   *   s \, f(R s) = s_0 .
   * \f]
   *
   * The left hand side is zero at \f$ s = 0 \f$ and strictly increasing, so the root is unique,
   * and since \f$ f \ge 1 \f$ it is bracketed by \f$ [0, s_0] \f$. This is solved by a Newton
   * iteration safeguarded by that bracket.
   *
   * Solving here rather than reading a drag functor is what keeps the functor dependency graph
   * acyclic: a drag material formed from the slip velocity, which is the correct definition,
   * cannot also be an input to the slip velocity.
   *
   * @param stokes_speed the slip magnitude in the Stokes limit, \f$ s_0 \f$
   * @param reynolds_per_speed the factor \f$ R \f$ converting a speed into a Reynolds number
   */
  static Real solveSlipSpeed(Real stokes_speed, Real reynolds_per_speed);

  /**
   * The factor converting the slip velocity into the diffusion velocity, \f$ 1 - c_d \f$.
   *
   * Written as \f$ (\rho_m - \alpha \rho_d) / \rho_m \f$ so that only the mixture density,
   * the dispersed phase density and the phase fraction are needed. The phase fraction is clamped
   * into [0, 1], matching the clamping applied by the mixture property material.
   */
  template <typename SpaceArg, typename StateArg>
  Real diffusionVelocityFactor(const SpaceArg & r, const StateArg & t) const
  {
    const auto rho_m = _rho_mixture(r, t);
    if (rho_m <= 0.0)
      return 1.0;
    const auto fd = std::clamp(_f_d(r, t), 0.0, 1.0);
    return (rho_m - fd * _rho_d(r, t)) / rho_m;
  }

  /// Dimension of the domain
  const unsigned int _dim;

  /// Velocity in the x direction
  MooseLinearVariableFVReal * const _u_var;
  /// Velocity in the y direction
  MooseLinearVariableFVReal * const _v_var;
  /// Velocity in the z direction
  MooseLinearVariableFVReal * const _w_var;

  /// Mixture density
  const Moose::Functor<Real> & _rho_mixture;
  /// Dispersed phase density
  const Moose::Functor<Real> & _rho_d;
  /// Continuous phase dynamic viscosity. Both the particle relaxation time and the particle
  /// Reynolds number are formed from it; see the note in validParams on why this is not the
  /// mixture viscosity
  const Moose::Functor<Real> & _mu_c;

  /// Volume fraction of the dispersed phase, needed to convert the slip velocity into the
  /// diffusion velocity
  const Moose::Functor<Real> & _f_d;

  /// Gravity acceleration vector
  const RealVectorValue _gravity;
  /// Body force scaling
  const Real _force_scale;
  /// Body force function
  const Function & _force_function;
  /// Body force postprocessor
  const PostprocessorValue & _force_postprocessor;
  /// Body force direction
  const RealVectorValue _force_direction;

  /// Prescribed linear drag function. Null when the drag is computed internally from the
  /// Schiller and Naumann correlation, see solveSlipSpeed
  const Moose::Functor<Real> * const _linear_friction;

  /// Continuous phase density, used to form the particle Reynolds number
  const Moose::Functor<Real> * const _rho_c;

  /// Which drag law to close the force balance with
  const DragModelEnum _drag_model;

  /// Surface tension, needed by the distorted particle drag
  const Moose::Functor<Real> * const _sigma;

  /// Exponent p of the hindrance factor (1 - alpha)^p applied to the slip velocity. Zero leaves
  /// the single-particle result untouched
  const Real _swarm_exponent;

  /// Frictional pressure gradients of the two phase and single particle systems, M_F and M_Finf,
  /// equations (41) and (24) of the Hibiki and Ishii paper cited on ishiiZuberSlipSpeed below.
  /// Only the 'ishii-zuber' drag model uses them.
  const Moose::Functor<Real> & _friction_pressure_gradient;
  const Moose::Functor<Real> & _single_particle_friction_pressure_gradient;

  /**
   * The relative velocity of the distorted particle multi-bubble system.
   *
   * The drag correlation this is named for is
   *
   *   Ishii, M. and Zuber, N., "Drag coefficient and relative velocity in bubbly, droplet or
   *   particulate flows", AIChE Journal 25 (1979) 843-855,
   *
   * and the form implemented here, which carries the frictional pressure gradients, is taken from
   *
   *   Hibiki, T. and Ishii, M., "One-dimensional drift-flux model and constitutive equations for
   *   relative motion between phases in various two-phase flow regimes", International Journal of
   *   Heat and Mass Transfer 46 (2003) 4935-4948.
   *
   * Every equation number quoted here and in the implementation refers to that second paper.
   * This function is its equations (45) to (50), which its equation (60) and Ishii's
   * (1 - alpha)^1.75 correlation are a power law fit to. Returns the magnitude of the slip, from
   * which the drift velocity follows as (1 - alpha) times it.
   *
   * @param alpha volume fraction of the dispersed phase
   * @param buoyancy the driving term, |rho_c - rho_d| times the acceleration, standing for the
   *        Delta rho g_z of the reference
   * @param m_f the frictional pressure gradient of the two phase flow
   * @param m_f_inf the frictional pressure gradient of the single particle system
   * @param sigma surface tension
   * @param rho_c density of the continuous phase
   */
  static Real ishiiZuberSlipSpeed(
      Real alpha, Real buoyancy, Real m_f, Real m_f_inf, Real sigma, Real rho_c);
  /// Particle diameter in the dispersed phase
  const Moose::Functor<Real> & _particle_diameter;

  /// Index of the velocity component x|y|z
  const unsigned int _index;
};
