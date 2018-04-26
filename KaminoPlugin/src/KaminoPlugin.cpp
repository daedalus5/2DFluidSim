#include <UT/UT_DSOVersion.h>
#include <UT/UT_Math.h>
#include <UT/UT_Interrupt.h>
#include <GU/GU_Detail.h>
#include <GU/GU_PrimPoly.h>
#include <CH/CH_LocalVariable.h>
#include <PRM/PRM_Include.h>
#include <PRM/PRM_SpareData.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>

#include <limits.h>
#include "KaminoPlugin.h"
using namespace HDK_Kamino;

///
/// newSopOperator is the hook that Houdini grabs from this dll
/// and invokes to register the SOP.  In this case we add ourselves
/// to the specified operator table.
///
void
newSopOperator(OP_OperatorTable *table)
{
    table->addOperator(
	    new OP_Operator("KaminoPlugin",			// Internal name
			    "Kamino",						// UI name
			     SOP_Kamino::myConstructor,		// How to build the SOP
			     SOP_Kamino::myTemplateList,	// My parameters
			     0,								// Min # of sources
			     0,								// Max # of sources
			     SOP_Kamino::myVariables,		// Local variables
			     OP_FLAG_GENERATOR)				// Flag it as generator
	    );
}

static enum Params {radius, nTheta, particleDensity, dt, DT, frames};

static PRM_Name
generateCommandName("generateCommand", "Generate");

static PRM_Name names[] =
{
	PRM_Name("radius", "Radius"),
	PRM_Name("nTheta", "Number of Latitudinal Subdivisions"),
	PRM_Name("particleDensity", "Max Particle Density"),
	PRM_Name("dt", "Time Step"),
	PRM_Name("DT", "Frame Rate"),
	PRM_Name("frames", "Number of Frames"),
};

static PRM_Default defaultParams[] =
{
	PRM_Default(5.0),
	PRM_Default(64),
	PRM_Default(100.0),
	PRM_Default(0.005),
	PRM_Default(0.041666667),
	PRM_Default(1000)
};

PRM_Template
SOP_Kamino::myTemplateList[] = 
{
	PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, names + radius, defaultParams + radius, 0),
	PRM_Template(PRM_INT, PRM_Template::PRM_EXPORT_MIN, 1, names + nTheta, defaultParams + nTheta, 0),
	PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, names + particleDensity, defaultParams + particleDensity, 0),
	PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, names + dt, defaultParams + dt, 0),
	PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, names + DT, defaultParams + DT, 0),
	PRM_Template(PRM_INT, PRM_Template::PRM_EXPORT_MIN, 1, names + frames, defaultParams + frames, 0),
	//PRM_Template(PRM_STRING, PRM_Template::PRM_EXPORT_MIN, 1, names + gram, defaultParams + gram, 0),
	PRM_Template(PRM_CALLBACK, 1, &generateCommandName, 0, 0, 0, Kamino::run),
    PRM_Template()
};


// Here's how we define local variables for the SOP.
enum {
	VAR_PT,		// Point number of the star
	VAR_NPT		// Number of points in the star
};

CH_LocalVariable
SOP_Kamino::myVariables[] = {
    { "PT", VAR_PT, 0 },		// The table provides a mapping
    { "NPT", VAR_NPT, 0 },		// from text string to integer token
    { 0, 0, 0 },
};

bool
SOP_Kamino::evalVariableValue(fpreal &val, int index, int thread)
{
    // myCurrPoint will be negative when we're not cooking so only try to
    // handle the local variables when we have a valid myCurrPoint index.
    if (myCurrPoint >= 0)
    {
	// Note that "gdp" may be null here, so we do the safe thing
	// and cache values we are interested in.
	switch (index)
	{
	    case VAR_PT:
		val = (fpreal) myCurrPoint;
		return true;
	    case VAR_NPT:
		val = (fpreal) myTotalPoints;
		return true;
	    default:
		/* do nothing */;
	}
    }
    // Not one of our variables, must delegate to the base class.
    return SOP_Node::evalVariableValue(val, index, thread);
}

OP_Node *
SOP_Kamino::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
    return new SOP_Kamino(net, name, op);
}

SOP_Kamino::SOP_Kamino(OP_Network *net, const char *name, OP_Operator *op)
	: SOP_Node(net, name, op)
{
    myCurrPoint = -1;	// To prevent garbage values from being returned
}

SOP_Kamino::~SOP_Kamino() {}

unsigned
SOP_Kamino::disableParms()
{
    return 0;
}

OP_ERROR
SOP_Kamino::cookMySop(OP_Context &context)
{
	fpreal now = context.getTime();

	fpreal rad = this->getRadius(now);
	exint ntheta = this->getNTheta(now);
	fpreal dens = this->getDensity(now);
	fpreal timestep = this->getdt(now);
	fpreal frameRate = this->getDT(now);
	exint frameCount = this->getFrames(now);

	// PUT YOUR CODE HERE
	int			 divisions;
	UT_Interrupt	*boss;

	// Since we don't have inputs, we don't need to lock them.

	divisions = 5;				// We need twice our divisions of points
	myTotalPoints = divisions;	// Set the NPT local variable value
	myCurrPoint = 0;			// Initialize the PT local variable
	
	// Check to see that there hasn't been a critical error in cooking the SOP.
	if (error() < UT_ERROR_ABORT)
	{
		boss = UTgetInterrupt();
		if (divisions < 4)
		{
			// With the range restriction we have on the divisions, this
			//	is actually impossible, but it shows how to add an error
			//	message or warning to the SOP.
			addWarning(SOP_MESSAGE, "Invalid divisions");
			divisions = 4;
		}
		gdp->clearAndDestroy();

		// Start the interrupt server
		if (boss->opStart("Running Solver"))
		{
			// PUT YOUR CODE HERE

			////////////////////////////////////////////////////////////////////////////////////////////

			select(GU_SPrimitive);
		}

		// Tell the interrupt server that we've completed. Must do this
		// regardless of what opStart() returns.
		boss->opEnd();
	}

	myCurrPoint = -1;
    return error();
}

