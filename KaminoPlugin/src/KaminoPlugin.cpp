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

#include <UT/UT_NTStreamUtil.h>
#include <UT/UT_IStream.h>
#include <CMD/CMD_Args.h>
#include <PI/PI_ResourceManager.h>
#include <MOT/MOT_Director.h>
#include <CMD/CMD_Manager.h>

#include <fstream>

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

/*
* Declare the GUI parameters to the nodes
*/

static enum Params {radius, nTheta, particleDensity, dt, DT, frames, densityImage, solidImage, colorImage, spin, uphi, utheta, vphi, vtheta};

static PRM_Name names[] =
{
	PRM_Name("radius", "Radius"),
	PRM_Name("nTheta", "Number of Latitudinal Subdivisions"),
	PRM_Name("particleDensity", "Max Particle Density"),
	PRM_Name("dt", "Time Step"),
	PRM_Name("DT", "Frame Rate"),
	PRM_Name("frames", "Number of Frames"),
	PRM_Name("spin", "Spin Speed"),
	PRM_Name("uphi", "u Variation in Longitude"),
	PRM_Name("utheta", "u Variation in latitude"),
	PRM_Name("vphi", "v Variation in Longitude"),
	PRM_Name("vtheta", "v Variation in Latitude"),
	PRM_Name("densityImage", "Density File"),
	PRM_Name("solidImage", "Solid File"),
	PRM_Name("colorImage", "Color Image File"),
};

static PRM_Default defaultParams[] =
{
	PRM_Default(5.0),
	PRM_Default(64),
	PRM_Default(100.0),
	PRM_Default(0.005),
	PRM_Default(0.041666667),
	PRM_Default(1000),
	PRM_Default(0.0),
	PRM_Default(1),
	PRM_Default(1),
	PRM_Default(1),
	PRM_Default(1),
	PRM_Default(0.0, ""),
	PRM_Default(0.0, ""),
	PRM_Default(0.0, ""),
};

static PRM_Name generateCommandName("generateCommand", "Run");

PRM_Template SOP_Kamino::myTemplateList[] = 
{
	PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, names + radius, defaultParams + radius, 0),
	PRM_Template(PRM_INT, PRM_Template::PRM_EXPORT_MIN, 1, names + nTheta, defaultParams + nTheta, 0),
	PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, names + particleDensity, defaultParams + particleDensity, 0),
	PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, names + dt, defaultParams + dt, 0),
	PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, names + DT, defaultParams + DT, 0),
	PRM_Template(PRM_INT, PRM_Template::PRM_EXPORT_MIN, 1, names + frames, defaultParams + frames, 0),
	PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, names + spin, defaultParams + spin, 0),
	PRM_Template(PRM_INT, PRM_Template::PRM_EXPORT_MIN, 1, names + uphi, defaultParams + uphi, 0),
	PRM_Template(PRM_INT, PRM_Template::PRM_EXPORT_MIN, 1, names + utheta, defaultParams + utheta, 0),
	PRM_Template(PRM_INT, PRM_Template::PRM_EXPORT_MIN, 1, names + vphi, defaultParams + vphi, 0),
	PRM_Template(PRM_INT, PRM_Template::PRM_EXPORT_MIN, 1, names + vtheta, defaultParams + vtheta, 0),
	PRM_Template(PRM_STRING, PRM_Template::PRM_EXPORT_MIN, 1, names + densityImage, defaultParams + densityImage, 0),
	PRM_Template(PRM_STRING, PRM_Template::PRM_EXPORT_MIN, 1, names + solidImage, defaultParams + solidImage, 0),
	PRM_Template(PRM_STRING, PRM_Template::PRM_EXPORT_MIN, 1, names + colorImage, defaultParams + colorImage, 0),
	PRM_Template(PRM_CALLBACK, 1, &generateCommandName, 0, 0, 0, SOP_Kamino::generateCallBack),
    PRM_Template()
};

std::string SOP_Kamino::pluginDir = std::getenv("CUSTOM_DSO_PATH");
static CMD_Manager m("Name");

int SOP_Kamino::generateCallBack(void* data, int index, float time, const PRM_Template*)
{
	SOP_Kamino* theNode = (SOP_Kamino*)data;
	std::string cmd = "system ";
	cmd = cmd + pluginDir + "/KaminoCoreSolver.exe ";
	cmd = cmd + pluginDir + "/configKamino.txt";
	UT_String hCmd = cmd.c_str();
	OPgetDirector()->getCommandManager()->execute(hCmd);
	
	m.execute(hCmd);
	return 1;
}

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
	: SOP_Node(net, name, op)//, pointer2Kamino(nullptr)
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
	fpreal A = this->getA(now);
	exint B = this->getB(now);
	exint C = this->getC(now);
	exint D = this->getD(now);
	exint E = this->getE(now);
	
	UT_String temp = "";
	this->getDensityFile(temp, now);
	std::string densityFile = temp.toStdString();

	temp = "";
	this->getSolidFile(temp, now);
	std::string solidFile = temp.toStdString();

	temp = "";
	this->getColorFile(temp, now);
	std::string colorFile = temp.toStdString();

	std::fstream fout;
	//std::string dsoPath = std::getenv("CUSTOM_DSO_PATH");
	fout.open(this->pluginDir + "/configKamino.txt", std::ios::out);
	fout << rad << std::endl;
	fout << ntheta << std::endl;
	fout << dens << std::endl;
	fout << timestep << std::endl;
	fout << frameRate << std::endl;
	fout << frameCount << std::endl;
	fout << A << std::endl;
	fout << B << std::endl;
	fout << C << std::endl;
	fout << D << std::endl;
	fout << E << std::endl;
	fout << this->pluginDir + "/output/frame" << std::endl;
	fout << this->pluginDir + "/particles/frame" << std::endl;
	if (densityFile.size() == 0)
	{
		fout << "null" << std::endl;
	}
	else
	{
		fout << densityFile << std::endl;
	}
	if (solidFile.size() == 0)
	{
		fout << "null" << std::endl;
	}
	else
	{
		fout << solidFile << std::endl;
	}
	if (colorFile.size() == 0)
	{
		fout << "null" << std::endl;
	}
	else
	{
		fout << colorFile << std::endl;
	}
	fout.close();

	OP_AutoLockInputs inputs(this);
	// Check if locking caused an error
	if (inputs.lock(context) >= UT_ERROR_ABORT)
		return error();

	OP_Node::flags().timeDep = 1;

    return error();
}

