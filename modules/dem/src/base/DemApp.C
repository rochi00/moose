#include "DemApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
DemApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

DemApp::DemApp(const InputParameters & parameters) : MooseApp(parameters)
{
  DemApp::registerAll(_factory, _action_factory, _syntax);
}

DemApp::~DemApp() {}

void
DemApp::registerAll(Factory & f, ActionFactory & af, Syntax & /*s*/)
{
  Registry::registerObjectsTo(f, {"DemApp"});
  Registry::registerActionsTo(af, {"DemApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
DemApp::registerApps()
{
  const std::string doc = "ArborX bounding volume hierarchy broad phase for the DEM neighbor list ";
#ifdef DEM_HAVE_ARBORX
  addBoolCapability("arborx", true, doc + "is available.");
#else
  addBoolCapability("arborx", false, doc + "is not available.");
#endif

  registerApp(DemApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
DemApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  DemApp::registerAll(f, af, s);
}
extern "C" void
DemApp__registerApps()
{
  DemApp::registerApps();
}
