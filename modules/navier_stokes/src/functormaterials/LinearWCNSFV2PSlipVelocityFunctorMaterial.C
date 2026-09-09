//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearWCNSFV2PSlipVelocityFunctorMaterial.h"
#include "FEProblemBase.h"
#include "Function.h"
#include "NS.h"
#include "NavierStokesMethods.h"

#include "libmesh/utility.h"

#include <limits>

registerMooseObject("NavierStokesApp", LinearWCNSFV2PSlipVelocityFunctorMaterial);

InputParameters
LinearWCNSFV2PSlipVelocityFunctorMaterial::validParams()
{
  InputParameters params = FunctorMaterial::validParams();
  params.addClassDescription(
      "Computes the slip velocity, and optionally the diffusion (drift) velocity derived from it, "
      "for the two-phase mixture model using the linear finite volume discretization.");
  params.addRequiredParam<SolverVariableName>("u", "The velocity in the x direction.");
  params.addParam<SolverVariableName>("v", "The velocity in the y direction.");
  params.addParam<SolverVariableName>("w", "The velocity in the z direction.");
  params.addRequiredParam<MooseFunctorName>(
      NS::density,
      "Mixture density. This is the density that appears in the buoyancy factor "
      "(rho_d - rho_m) / rho_d of the algebraic slip closure, see VTT Publications 288 equation "
      "(58), not the continuous phase density.");
  params.addRequiredParam<MooseFunctorName>("rho_d", "Dispersed phase density.");
  params.addRequiredParam<MooseFunctorName>(
      NS::mu,
      "Continuous phase dynamic viscosity. Both the particle relaxation time and the particle "
      "Reynolds number are formed from it, which is the single particle form of the closure. This "
      "is deliberately not the mixture viscosity. Manninen substitutes an apparent mixture "
      "viscosity here as a way of folding the swarm correction into the drag, but that "
      "substitution is only valid with the apparent viscosity of Ishii and Zuber, VTT "
      "Publications 288 equation (43), which rises with the dispersed fraction; a volume "
      "averaged mixture viscosity falls with it and would invert the correction. Concentration "
      "effects belong in 'swarm_exponent' or in the 'ishii-zuber' drag model instead.");
  params.addParam<RealVectorValue>(
      "gravity", RealVectorValue(0, 0, 0), "Gravity acceleration vector");
  params.addParam<Real>("force_value", 0.0, "Coefficient to multiply by the body force term");
  params.addParam<FunctionName>("force_function", "0", "A function that describes the body force");
  params.addParam<PostprocessorName>(
      "force_postprocessor", 0, "A postprocessor whose value is multiplied by the body force");
  params.addParam<RealVectorValue>(
      "force_direction", RealVectorValue(1, 0, 0), "Direction of the body force");
  params.addParam<MooseFunctorName>(
      "linear_coef_name",
      "Prescribed linear drag function. Supply this to impose the drag directly. Leave it unset "
      "and set 'use_dispersed_phase_drag_model' to have the Schiller and Naumann correlation "
      "solved together with the force balance instead.");
  params.addParam<bool>(
      "use_dispersed_phase_drag_model",
      false,
      "Whether to evaluate the drag from the Schiller and Naumann correlation, solved "
      "self-consistently with the slip velocity so that the particle Reynolds number is formed "
      "from the slip velocity as its definition requires.");
  params.addParam<MooseFunctorName>(
      "rho_c", "Continuous phase density, used to form the particle Reynolds number.");
  MooseEnum drag_model("schiller-naumann distorted-particle automatic ishii-zuber",
                       "schiller-naumann");
  params.addParam<MooseEnum>(
      "drag_model",
      drag_model,
      "Drag law closing the force balance when 'use_dispersed_phase_drag_model' is set. "
      "'schiller-naumann' is the rigid sphere law and is accurate while the particle stays "
      "spherical. 'distorted-particle' is the deformed bubble law, whose terminal velocity is "
      "independent of the particle size and reproduces Ishii's drift velocity correlation. "
      "'automatic' takes whichever of the two resists more, which selects the regime by itself. "
      "'ishii-zuber' is the multi-bubble relation of Ishii and Zuber, which predicts the relative "
      "velocity directly instead of supplying a drag to close the balance; it carries its own "
      "concentration dependence and therefore ignores 'swarm_exponent'.");
  params.addParam<MooseFunctorName>(
      "friction_pressure_gradient",
      "0",
      "The frictional pressure gradient of the two phase flow, M_F = 4 tau_fw / D = (-dp/dz)_F, "
      "which is equation (41) of Hibiki and Ishii, \"One-dimensional drift-flux model and "
      "constitutive equations for relative motion between phases in various two-phase flow "
      "regimes\", Int. J. Heat Mass Transfer 46 (2003) 4935-4948. Only used by the 'ishii-zuber' "
      "drag model, where it carries the effect of the wall friction on the relative velocity. "
      "Zero, the default, leaves the balance gravity dominant, which is the limit in which that "
      "model reduces to Ishii's correlation. That paper obtains it from the correlation of "
      "Lockhart and Martinelli.");
  params.addParam<MooseFunctorName>(
      "single_particle_friction_pressure_gradient",
      "0",
      "The frictional pressure gradient of the corresponding single particle system, "
      "M_Finf = (f / 2D) rho_f <v_finf>^2. This is equation (24) of the Hibiki and Ishii paper "
      "cited under 'friction_pressure_gradient'; it enters both the terminal velocity of that "
      "paper's equation (49) and the ratio of its equation (46), so it is not simply a correction "
      "on the hindrance.");
  params.addParam<MooseFunctorName>(
      "surface_tension",
      "Surface tension between the phases. Required by the distorted particle drag.");
  params.addParam<Real>(
      "swarm_exponent",
      0.0,
      "Exponent p of the hindrance factor (1 - alpha)^p multiplying the slip velocity. A particle "
      "in a swarm meets more resistance than an isolated one, and 'schiller-naumann' and "
      "'distorted-particle' are both single particle laws, so they need this correction. Zero, "
      "the default, leaves the single particle result untouched. A value of 0.75 makes the "
      "void-fraction-weighted mean drift velocity scale as (1 - alpha)^1.75, which is the swarm "
      "dependence of Ishii's correlation, since the remaining factor of (1 - alpha) is kinematic. "
      "It is ignored by 'ishii-zuber', which is a multi-particle correlation and carries its own "
      "concentration dependence.");

  params.addParam<MooseFunctorName>(
      "particle_diameter", 1.0, "Diameter of particles in the dispersed phase.");
  MooseEnum momentum_component("x=0 y=1 z=2");
  params.addRequiredParam<MooseEnum>(
      "momentum_component",
      momentum_component,
      "The component of the momentum equation that this material applies to.");
  params.addRequiredParam<MooseFunctorName>("slip_velocity_name", "the name of the slip velocity");
  params.addParam<MooseFunctorName>(
      "drift_velocity_name",
      "Name of the diffusion (drift) velocity property to declare, the velocity of the dispersed "
      "phase relative to the centre of mass of the mixture. Declared only when supplied.");
  params.addParam<MooseFunctorName>(
      "fd",
      0.0,
      "Volume fraction of the dispersed phase. It converts the slip velocity into the drift "
      "velocity, it sets the swarm hindrance factor, and it supplies the whole concentration "
      "dependence of the 'ishii-zuber' drag model.");
  params.renameParam("fd", "fraction_dispersed", "");
  params.addParam<unsigned short>("ghost_layers",
                                  3,
                                  "The number of layers of elements to ghost. With Rhie-Chow and "
                                  "the velocity gradient calculation below, we need 3");
  params.addRelationshipManager(
      "ElementSideNeighborLayers",
      Moose::RelationshipManagerType::GEOMETRIC | Moose::RelationshipManagerType::ALGEBRAIC |
          Moose::RelationshipManagerType::COUPLING,
      [](const InputParameters & obj_params, InputParameters & rm_params)
      {
        rm_params.set<unsigned short>("layers") = obj_params.get<unsigned short>("ghost_layers");
      });
  return params;
}

MooseLinearVariableFVReal &
LinearWCNSFV2PSlipVelocityFunctorMaterial::getVelocityVariable(const std::string & param_name)
{
  auto * const var = dynamic_cast<MooseLinearVariableFVReal *>(
      &_fe_problem.getVariable(_tid, getParam<SolverVariableName>(param_name)));
  if (!var)
    paramError(param_name, "The velocity must be a MooseLinearVariableFVReal.");
  return *var;
}

LinearWCNSFV2PSlipVelocityFunctorMaterial::LinearWCNSFV2PSlipVelocityFunctorMaterial(
    const InputParameters & params)
  : FunctorMaterial(params),
    NonADFunctorInterface(this),
    _dim(_subproblem.mesh().dimension()),
    _u_var(&getVelocityVariable("u")),
    _v_var(_dim > 1 && isParamValid("v") ? &getVelocityVariable("v") : nullptr),
    _w_var(_dim > 2 && isParamValid("w") ? &getVelocityVariable("w") : nullptr),
    _rho_mixture(NonADFunctorInterface::getFunctor<Real>(NS::density)),
    _rho_d(NonADFunctorInterface::getFunctor<Real>("rho_d")),
    _mu_c(NonADFunctorInterface::getFunctor<Real>(NS::mu)),
    _f_d(NonADFunctorInterface::getFunctor<Real>("fd")),
    _gravity(getParam<RealVectorValue>("gravity")),
    _force_scale(getParam<Real>("force_value")),
    _force_function(getFunction("force_function")),
    _force_postprocessor(getPostprocessorValue("force_postprocessor")),
    _force_direction(getParam<RealVectorValue>("force_direction")),
    _linear_friction(isParamValid("linear_coef_name")
                         ? &NonADFunctorInterface::getFunctor<Real>("linear_coef_name")
                         : nullptr),
    _rho_c(isParamValid("rho_c") ? &NonADFunctorInterface::getFunctor<Real>("rho_c") : nullptr),
    _drag_model(getParam<MooseEnum>("drag_model").getEnum<DragModelEnum>()),
    _sigma(isParamValid("surface_tension")
               ? &NonADFunctorInterface::getFunctor<Real>("surface_tension")
               : nullptr),
    _swarm_exponent(getParam<Real>("swarm_exponent")),
    _friction_pressure_gradient(
        NonADFunctorInterface::getFunctor<Real>("friction_pressure_gradient")),
    _single_particle_friction_pressure_gradient(
        NonADFunctorInterface::getFunctor<Real>("single_particle_friction_pressure_gradient")),
    _particle_diameter(NonADFunctorInterface::getFunctor<Real>("particle_diameter")),
    _index(getParam<MooseEnum>("momentum_component"))
{
  const bool internal_drag = getParam<bool>("use_dispersed_phase_drag_model");
  if (internal_drag && _linear_friction)
    paramError("linear_coef_name",
               "A prescribed drag function cannot be combined with "
               "'use_dispersed_phase_drag_model', which computes the drag internally.");
  if (!internal_drag && !_linear_friction)
    paramError("linear_coef_name",
               "Either supply a drag function here, or set 'use_dispersed_phase_drag_model' to "
               "compute it from the Schiller and Naumann correlation.");
  if (internal_drag && !_rho_c)
    paramError("use_dispersed_phase_drag_model",
               "The continuous phase density 'rho_c' is required to form the particle Reynolds "
               "number of the drag correlation.");
  if (internal_drag && _drag_model != DragModelEnum::SCHILLER_NAUMANN)
  {
    if (!_sigma)
      paramError("surface_tension",
                 "The 'distorted-particle' and 'ishii-zuber' drag models are both set by the "
                 "balance of buoyancy against surface tension, so the surface tension is "
                 "required.");
    // Only the distorted particle law is scaled by the gravity magnitude itself. The Ishii and
    // Zuber relation drives its buoyancy from the full particle acceleration, so a case driven by
    // a body force rather than by gravity is legitimate for it.
    if (_gravity.norm() == 0 && _drag_model != DragModelEnum::ISHII_ZUBER)
      paramError("gravity",
                 "The 'distorted-particle' drag is scaled by g * delta_rho / sigma and is not "
                 "defined without gravity.");
  }

  // The phase fraction defaults to zero, which silently disables both of the mechanisms that
  // depend on it rather than failing, so require it wherever it actually carries the physics.
  if (!isParamSetByUser("fraction_dispersed"))
  {
    if (_swarm_exponent != 0.0)
      paramError("fraction_dispersed",
                 "A swarm exponent was supplied but the dispersed phase fraction was not. The "
                 "hindrance factor is (1 - alpha)^p, so without the phase fraction it is exactly "
                 "one and the swarm correction requested here would have no effect.");
    if (_drag_model == DragModelEnum::ISHII_ZUBER)
      paramError("fraction_dispersed",
                 "The 'ishii-zuber' drag model is a multi-particle correlation whose entire "
                 "concentration dependence comes from the dispersed phase fraction, so it must be "
                 "supplied. Without it the model returns the isolated particle result.");
  }

  if (_dim > 1 && !isParamValid("v"))
    paramError("v", "In two or more dimensions, the v velocity must be supplied.");
  if (_dim > 2 && !isParamValid("w"))
    paramError("w", "In three dimensions, the w velocity must be supplied.");

  // The advective part of the particle acceleration needs the velocity gradients. Requesting them
  // here is unconditional, unlike the nonlinear counterpart which only reaches this call when a
  // cast succeeds and silently skips it otherwise.
  _u_var->computeCellGradients();
  if (_v_var)
    _v_var->computeCellGradients();
  if (_w_var)
    _w_var->computeCellGradients();

  const auto & slip_velocity = addFunctorProperty<Real>(
      getParam<MooseFunctorName>("slip_velocity_name"),
      [this](const auto & r, const auto & t) -> Real
      {
        // Guards the division below against a prescribed friction function evaluating to zero
        constexpr Real offset = 1e-15;

        RealVectorValue term_advection(0, 0, 0);
        RealVectorValue term_transient(0, 0, 0);
        const RealVectorValue term_force(
            _force_scale * _force_postprocessor *
            _force_function.value(_t, _current_elem->vertex_average()) * _force_direction);

        // The time derivative is taken from the variable itself, which routes it through the
        // problem's time integrator, see MooseLinearVariableFV::evaluateDot. That matters because
        // this acceleration feeds a slip velocity which in turn feeds the momentum equation, so a
        // hand rolled first order difference here would drop the whole scheme to first order
        // under, for instance, BDF2.
        if (_subproblem.isTransient())
        {
          term_transient(0) = raw_value(_u_var->dot(r, t));
          if (_v_var)
            term_transient(1) = raw_value(_v_var->dot(r, t));
          if (_w_var)
            term_transient(2) = raw_value(_w_var->dot(r, t));
        }

        // Advective part of the material derivative, u . grad(u)
        const auto u_velocity = raw_value((*_u_var)(r, t));
        const auto u_grad = raw_value(_u_var->gradient(r, t));
        term_advection(0) += u_velocity * u_grad(0);
        if (_v_var)
        {
          const auto v_velocity = raw_value((*_v_var)(r, t));
          const auto v_grad = raw_value(_v_var->gradient(r, t));
          term_advection(0) += v_velocity * u_grad(1);
          term_advection(1) += u_velocity * v_grad(0) + v_velocity * v_grad(1);
          if (_w_var)
          {
            const auto w_velocity = raw_value((*_w_var)(r, t));
            const auto w_grad = raw_value(_w_var->gradient(r, t));
            term_advection(0) += w_velocity * u_grad(2);
            term_advection(1) += w_velocity * v_grad(2);
            term_advection(2) +=
                u_velocity * w_grad(0) + v_velocity * w_grad(1) + w_velocity * w_grad(2);
          }
        }

        // Manninen's algebraic slip closure, u_slip = tau_d / f_drag * (rho_d - rho_m)/rho_d * a,
        // with the particle relaxation time of Bilicki and Kestin.
        //
        // The buoyancy factor carries the MIXTURE density, because the pressure gradient the
        // particle sits in is the mixture's, while the relaxation time carries the CONTINUOUS
        // phase viscosity, because Stokes drag is exerted by the fluid the particle moves
        // through. That pairing is the single particle form of the closure.
        //
        // Manninen instead substitutes an apparent mixture viscosity into the relaxation time,
        // which is a way of folding the swarm correction into the drag itself. That substitution
        // is only valid with the Ishii and Zuber apparent viscosity of VTT Publications 288
        // equation (43), which rises with the dispersed fraction. A volume averaged mixture
        // viscosity falls with it and would invert the correction. This object therefore keeps
        // the drag single particle and leaves every concentration effect to swarm_exponent and to
        // the ishii-zuber drag model, so that the two cannot double count
        const Real density_scaling = (_rho_d(r, t) - _rho_mixture(r, t)) / _rho_d(r, t);
        const RealVectorValue acceleration_vec =
            -term_transient - term_advection + _gravity + term_force;

        const Real relaxation_time =
            _rho_d(r, t) * Utility::pow<2>(_particle_diameter(r, t)) / (18.0 * _mu_c(r, t));

        // The slip in the Stokes limit, f_drag = 1, which is also the whole answer when the drag
        // function is prescribed rather than computed
        const Real stokes_prefactor = relaxation_time * density_scaling;

        if (_linear_friction)
          return stokes_prefactor / ((*_linear_friction)(r, t) + offset) *
                 acceleration_vec(_index);

        // Otherwise solve the force balance and the drag correlation together, so that the
        // particle Reynolds number is formed from the slip velocity as its definition requires
        const Real stokes_speed = std::abs(stokes_prefactor) * acceleration_vec.norm();
        const Real rho_c = (*_rho_c)(r, t);
        const Real mu_c = _mu_c(r, t);

        // The Ishii and Zuber multi-bubble relation is a correlation for the relative velocity
        // itself, not a drag law, so it replaces the force balance rather than closing it. It also
        // carries its own concentration dependence, which is what the swarm exponent stands in for
        // elsewhere, so no hindrance is applied on top of it.
        if (_drag_model == DragModelEnum::ISHII_ZUBER)
        {
          const Real alpha = std::clamp(_f_d(r, t), 0.0, 1.0);
          const Real buoyancy = std::abs(rho_c - _rho_d(r, t)) * acceleration_vec.norm();
          const Real speed = ishiiZuberSlipSpeed(alpha,
                                                 buoyancy,
                                                 _friction_pressure_gradient(r, t),
                                                 _single_particle_friction_pressure_gradient(r, t),
                                                 (*_sigma)(r, t),
                                                 rho_c);
          // The correlation returns a magnitude. Its direction is that of the buoyancy, which is
          // the direction of the acceleration carrying the sign of rho_d - rho_m: a light dispersed
          // phase moves against the acceleration, a heavy one with it.
          const Real norm = acceleration_vec.norm();
          if (norm <= 0.0)
            return 0.0;
          return speed * std::copysign(1.0, density_scaling) * acceleration_vec(_index) / norm;
        }

        // Rigid sphere branch: s f(R s) = s0 with f monotone and at least unity
        Real slip_speed = std::numeric_limits<Real>::max();
        if (_drag_model != DragModelEnum::DISTORTED_PARTICLE)
        {
          const Real reynolds_per_speed =
              NS::particleReynoldsNumber(rho_c, _particle_diameter(r, t), 1.0, mu_c);
          slip_speed = solveSlipSpeed(stokes_speed, reynolds_per_speed);
        }

        // Distorted particle branch. Its drag function is linear in the slip speed, so the force
        // balance k s^2 = s0 has the closed form s = sqrt(s0 / k) and needs no iteration.
        if (_drag_model != DragModelEnum::SCHILLER_NAUMANN)
        {
          const Real k = NS::distortedDragFunctionPerSpeed(_particle_diameter(r, t),
                                                           rho_c,
                                                           mu_c,
                                                           std::abs(rho_c - _rho_d(r, t)),
                                                           (*_sigma)(r, t),
                                                           _gravity.norm());
          const Real distorted_speed = (k > 0.0) ? std::sqrt(stokes_speed / k) : stokes_speed;
          // A larger drag function gives a smaller slip, so taking the more resistant of the two
          // laws is the same as taking the smaller of the two speeds
          slip_speed = std::min(slip_speed, distorted_speed);
        }

        // Recover the drag function the chosen speed implies, and apply the closure with it. This
        // keeps a single expression for the returned component whichever branch was taken.
        const Real f_drag = (slip_speed > 0.0) ? stokes_speed / slip_speed : 1.0;

        // Hindrance of the swarm. Applied to the speed rather than to the drag function because
        // the two branches respond differently to a factor on the drag, being respectively linear
        // and square-root in it, whereas a factor on the speed means the same thing for both.
        const Real hindrance =
            (_swarm_exponent == 0.0)
                ? 1.0
                : std::pow(std::max(1.0 - std::clamp(_f_d(r, t), 0.0, 1.0), 0.0),
                           _swarm_exponent);

        return stokes_prefactor / f_drag * hindrance * acceleration_vec(_index);
      });

  // The diffusion (drift) velocity, the velocity of the dispersed phase relative to the centre of
  // mass of the mixture. This, not the slip velocity, is the velocity that appears in the phase
  // conservation equations, see VTT Publications 288 equation (28).
  if (isParamValid("drift_velocity_name"))
    addFunctorProperty<Real>(getParam<MooseFunctorName>("drift_velocity_name"),
                             [this, &slip_velocity](const auto & r, const auto & t) -> Real
                             { return diffusionVelocityFactor(r, t) * slip_velocity(r, t); });
}


Real
LinearWCNSFV2PSlipVelocityFunctorMaterial::solveSlipSpeed(const Real stokes_speed,
                                                          const Real reynolds_per_speed)
{
  // f(0) = 1, so a vanishing acceleration gives a vanishing slip and the drag never enters
  if (stokes_speed <= 0.0 || reynolds_per_speed <= 0.0)
    return stokes_speed;

  // Solve in Reynolds number rather than in speed. Multiplying s f(R s) = s0 through by R turns it
  // into Re f(Re) = B, with B = R s0 a dimensionless group formed entirely from inputs. The root is
  // unchanged; what is gained is that the branch of f can be chosen before iterating, from a
  // quantity that does not move as the iteration proceeds. The loop below therefore only ever
  // evaluates the Schiller and Naumann branch, and cannot step across the 0.2% seam dragFunction
  // takes at Re = 1000 the way an iteration re-testing its own iterate could.
  const Real driving_group = reynolds_per_speed * stokes_speed;

  // Below this the drag correction 0.15 Re^0.687 is 3e-13, smaller than the relative tolerance the
  // loop would converge to, so the Stokes answer is already the converged one.
  constexpr Real negligible_driving_group = 1e-17;
  if (driving_group < negligible_driving_group)
    return stokes_speed;

  // In Newton's regime f = 0.0183 Re, so the balance becomes 0.0183 Re^2 = B and is exact. The
  // branches of dragFunction change over at Re = 1000, which is this value of B.
  constexpr Real newton_regime_driving_group = 0.0183 * 1000.0 * 1000.0;
  if (driving_group >= newton_regime_driving_group)
    return std::sqrt(driving_group / 0.0183) / reynolds_per_speed;

  // g(Re) = Re f(Re) is zero at the origin and strictly increasing, and f >= 1 puts the root in
  // [0, B]. Newton is safeguarded by that bracket so that it cannot leave it.
  Real lower = 0.0;
  Real upper = driving_group;

  // One Picard step off the Stokes guess Re = B. It is exact in the Stokes limit and loosens
  // towards the transition, where the bracket absorbs the resulting overshoot in one pass.
  Real reynolds = driving_group / NS::schillerNaumannDragFunction(driving_group);

  // The residual falls below this relative tolerance in a handful of iterations; the cap is a
  // backstop, not the expected exit
  constexpr Real rel_tol = 1e-12;
  constexpr unsigned int max_its = 50;

  for ([[maybe_unused]] const auto it : make_range(max_its))
  {
    const Real drag = NS::schillerNaumannDragFunction(reynolds);
    const Real residual = reynolds * drag - driving_group;

    if (std::abs(residual) <= rel_tol * driving_group)
      break;

    if (residual > 0.0)
      upper = reynolds;
    else
      lower = reynolds;

    // d/dRe [Re f(Re)] = f(Re) + Re f'(Re), the second term written as the finite product
    const Real derivative = drag + NS::schillerNaumannDragDerivative(reynolds);
    const Real candidate = reynolds - residual / derivative;

    // Fall back on bisection if Newton steps outside the bracket
    reynolds = (candidate > lower && candidate < upper) ? candidate : 0.5 * (lower + upper);
  }

  return reynolds / reynolds_per_speed;
}


Real
LinearWCNSFV2PSlipVelocityFunctorMaterial::ishiiZuberSlipSpeed(const Real alpha,
                                                               const Real buoyancy,
                                                               const Real m_f,
                                                               const Real m_f_inf,
                                                               const Real sigma,
                                                               const Real rho_c)
{
  const Real one_minus = std::max(1.0 - alpha, 0.0);
  const Real denominator = buoyancy + m_f_inf;
  if (denominator <= 0.0 || sigma <= 0.0 || rho_c <= 0.0 || one_minus <= 0.0)
    return 0.0;

  // Equation numbers below are those of Hibiki and Ishii, "One-dimensional drift-flux model and
  // constitutive equations for relative motion between phases in various two-phase flow regimes",
  // Int. J. Heat Mass Transfer 46 (2003) 4935-4948, the paper this form is taken from. The
  // underlying drag correlation is that of Ishii and Zuber, AIChE Journal 25 (1979) 843-855.
  //
  // The terminal velocity of an isolated distorted particle, equation (49). The frictional
  // pressure gradient of the single particle system adds to the buoyancy here, so it raises the
  // terminal velocity rather than merely hindering it.
  const Real v_r_inf =
      std::sqrt(2.0) * std::pow(denominator * sigma / Utility::pow<2>(rho_c), 0.25);

  // The ratio of equation (46), through which the frictional pressure gradient of the two phase
  // flow enters. With both gradients zero it is simply 1 - alpha.
  const Real ratio = std::max((buoyancy * one_minus + m_f) / denominator, 0.0);

  // Equation (46) with the viscosity ratio of equation (47), mu_f / mu_m = 1 - alpha
  const Real f_alpha = one_minus * std::sqrt(ratio);

  // Equation (45), the multi-particle relative velocity. The constants are those of the Newton
  // regime drag ratio of Ishii and Zuber (1979), arranged so that the concentration factor is
  // unity at alpha = 0.
  const Real concentration =
      18.67 * f_alpha / (1.0 + 17.67 * std::pow(f_alpha, 6.0 / 7.0));

  // Equation (45) gives the drift velocity; the relative velocity carries one power of
  // (1 - alpha) less
  return v_r_inf * std::sqrt(one_minus) * concentration;
}
