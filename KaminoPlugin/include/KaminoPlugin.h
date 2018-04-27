#ifndef __LSYSTEM_PLUGIN_h__
#define __LSYSTEM_PLUGIN_h__

#include <SOP/SOP_Node.h>
#include "Kamino.h"

namespace HDK_Kamino {
class SOP_Kamino : public SOP_Node
{
public:
    static OP_Node *myConstructor(OP_Network*, const char *, OP_Operator *);

    /// Stores the description of the interface of the SOP in Houdini.
    /// Each parm template refers to a parameter.
    static PRM_Template myTemplateList[];

    /// This optional data stores the list of local variables.
    static CH_LocalVariable myVariables[];

protected:

	SOP_Kamino(OP_Network *net, const char *name, OP_Operator *op);
    virtual ~SOP_Kamino();

    /// Disable parameters according to other parameters.
    virtual unsigned disableParms();

    /// cookMySop does the actual work of the SOP computing
    virtual OP_ERROR cookMySop(OP_Context &context);

    /// This function is used to lookup local variables that you have
    /// defined specific to your SOP.
    virtual bool evalVariableValue(fpreal &val, int index, int thread);

    // Add virtual overload that delegates to the super class to avoid
    // shadow warnings.
    virtual bool evalVariableValue(UT_String &v, int i, int thread){
		return evalVariableValue(v, i, thread);
	}

private:
    /// The following list of accessors simplify evaluating the parameters
    /// of the SOP.
	fpreal getRadius(fpreal t)
	{
		return evalFloat("radius", 0, t);
	}
	exint getNTheta(fpreal t)
	{	
		return evalInt("nTheta", 0, t);
	}
	fpreal getDensity(fpreal t)
	{
		return evalFloat("particleDensity", 0, t);
	}
	fpreal getdt(fpreal t)
	{
		return evalFloat("dt", 0, t);
	}
	fpreal getDT(fpreal t)
	{
		return evalFloat("DT", 0, t);
	}
	exint getFrames(fpreal t)
	{
		return evalInt("frames", 0, t);
	}
	void getDensityFile(UT_String& str, fpreal t)
	{
		evalString(str, "densityImage", 0, t);
	}
	void getSolidFile(UT_String& str, fpreal t)
	{
		evalString(str, "solidImage", 0, t);
	}
	void getColorFile(UT_String& str, fpreal t)
	{
		evalString(str, "colorImage", 0, t);
	}


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Member variables are stored in the actual SOP, not with the geometry
    /// In this case these are just used to transfer data to the local 
    /// variable callback.
    /// Another use for local data is a cache to store expensive calculations.

	// NOTE : You can declare local variables here like this  
    int myCurrPoint;
    int myTotalPoints;
};
} // End HDK_Sample namespace

#endif
