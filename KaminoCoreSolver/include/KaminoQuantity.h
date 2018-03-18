# pragma once

# include <string>
# include <map>
# include <iostream>
# include <vector>
# ifndef _MSC_VER
# include "Partio.h"
# endif
# include <Eigen/Core>
# include <Eigen/Dense>
# include <cmath>

# define M_PI           3.14159265358979323846  /* pi */

# define DEBUGBUILD

// The solution to switch between double and float
typedef double fReal;


// Handy Lerp.
template <class Type>
Type KaminoLerp(const Type &fromEndPoint, const Type &toEndPoint, double factor);

// The attribute base class.
class KaminoQuantity
{
private:

	std::string attrName;

	/* Grid dimensions */
	/* These are only assigned by the Grid and not mutable */
	size_t nx;
	size_t ny;

	/* Grid size */
	fReal gridLen;

	/* Is this staggered? */
	fReal xOffset;
	fReal yOffset;

	/* Double buffer */
	fReal* thisStep;
	fReal* nextStep;

	/* Wrap things up */
	size_t getWarpedXIndex(fReal x);
	size_t getWarpedYIndex(fReal y);

	/* Get index */
	size_t getIndex(size_t x, size_t y);

public:
	/* Constructor */
	KaminoQuantity(std::string attributeName, size_t nx, size_t ny, fReal gridLen, fReal xOffset, fReal yOffset);
	/* Destructor */
	~KaminoQuantity();

	/* Swap the buffer */
	void swapBuffer();

	/* Get the current step */
	fReal getValueAt(size_t x, size_t y);
	/* Set the current step */
	void setValueAt(size_t x, size_t y, fReal val);
	/* Write to the next step */
	void writeValueTo(size_t x, size_t y, fReal val);
	/* Access */
	fReal& accessValueAt(size_t x, size_t y);
	/* Lerped Sampler using world coordinates */
	fReal sampleAt(fReal x, fReal y);
};


// The solver class.
class KaminoSolver
{
private:
	/* Grid dimensions */
	size_t nx;
	size_t ny;
	/* Grid size */
	fReal gridLen;

	/* So that it remembers all these attributes within */
	std::map<std::string, KaminoQuantity*> attributeTable;

	/* Something about time steps */
	fReal frameDuration;
	fReal timeStep;
	fReal timeElapsed;

	// TODO
	void advection();
	void projection();
	void bodyForce();

	// Swap all these buffers of the attributes.
	void swapAttrBuffers();

	/* distribute initial velocity values at grid points */
	void initialize_velocity();
	/* FBM noise function for velocity distribution */
	fReal FBM(const fReal x, const fReal y);
	/* 2D noise interpolation function for smooth FBM noise */
	fReal interpNoise2D(const fReal x, const fReal y) const;
	/* returns a pseudorandom number between -1 and 1 */
	fReal rand(const Eigen::Matrix<fReal, 2, 1> vecA) const;

public:
	KaminoSolver(size_t nx, size_t ny, fReal gridLength, fReal frameDuration = 1.0 / 30.0);
	~KaminoSolver();

	void stepForward(fReal timeStep);

	void addAttr(std::string name, fReal xOffset = 0.0, fReal yOffset = 0.0);
	
	KaminoQuantity* getAttributeNamed(std::string name);
	KaminoQuantity* operator[](std::string name);

	void write_data_bgeo(const std::string& s, const int frame);
};