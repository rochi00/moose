//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "DemTestApp.h"
#include "DemApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
DemTestApp::validParams()
{
  InputParameters params = DemApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

DemTestApp::DemTestApp(const InputParameters & parameters) : MooseApp(parameters)
{
  DemTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

DemTestApp::~DemTestApp() {}

void
DemTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  DemApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"DemTestApp"});
    Registry::registerActionsTo(af, {"DemTestApp"});
  }
}

void
DemTestApp::registerApps()
{
  registerApp(DemApp);
  registerApp(DemTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
DemTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  DemTestApp::registerAll(f, af, s);
}
extern "C" void
DemTestApp__registerApps()
{
  DemTestApp::registerApps();
}
