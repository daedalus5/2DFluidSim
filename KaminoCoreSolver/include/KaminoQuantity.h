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
# include <Eigen/IterativeLinearSolvers>
# include <unsupported/Eigen/IterativeSolvers>
# include <boost/math/tools/roots.hpp>

# define M_PI           3.14159265358979323846  /* pi */

# define DEBUGBUILD

// The solution to switch between double and float
typedef double fReal;


// Handy Lerp.
template <class Type>
Type KaminoLerp(const Type &fromEndPoint, const Type &toEndPoint, double factor)
{
	return (1.0 - factor) * fromEndPoint + factor * toEndPoint;
}

enum gridType { FLUIDGRID, SOLIDGRID, AIRGRID };

// The attribute base class.
class KaminoQuantity
{
private:

	std::string attrName;

	/* Grid dimensions */
	/* These are only assigned by the Grid and not mutable */
	size_t nPhi;
	size_t nTheta;

	/* Grid size */
	fReal gridLen;
	fReal invGridLen;

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
	inline size_t getIndex(size_t x, size_t y);

public:
	/* Constructor */
	KaminoQuantity(std::string attributeName, size_t nx, size_t ny, fReal gridLen, fReal xOffset, fReal yOffset);
	/* Destructor */
	~KaminoQuantity();

	/* Swap the buffer */
	void swapBuffer();

	/* Get nx */
	size_t getNPhi();
	/* Get ny */
	size_t getNTheta();
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
	/* Given the index, show its origin in world coordinates*/
	fReal getXCoordAtIndex(size_t x);
	fReal getYCoordAtIndex(size_t y);
};

// The solver class.
class KaminoSolver
{
private:
	/* Grid types */
	gridType* gridTypes;
	/* Grid dimensions */
	size_t nPhi;
	size_t nTheta;
	/* Radius of sphere */
	fReal radius;
	/* Grid size */
	fReal gridLen;
	/* Laplacian Matrix */
	Eigen::SparseMatrix<fReal> Laplacian;

	/* So that it remembers all these attributes within */
	std::map<std::string, KaminoQuantity*> attributeTable;

	/* Something about time steps */
	fReal frameDuration;
	fReal timeStep;
	fReal timeElapsed;

	gridType getGridTypeAt(size_t x, size_t y);

	void advection();
	void geometric();
	void projection();
	void bodyForce();

	// Swap all these buffers of the attributes.
	void swapAttrBuffers();

	/* distribute initial velocity values at grid points */
	void initialize_velocity();
	/* initialize pressure attribute */
	void initialize_pressure();
	/* initialize test case */
	void initialize_test();
	/* which grids are solid? */
	void initialize_boundary();
	/* FBM noise function for velocity distribution */
	fReal FBM(const fReal x, const fReal y);
	/* 2D noise interpolation function for smooth FBM noise */
	fReal interpNoise2D(const fReal x, const fReal y) const;
	/* returns a pseudorandom number between -1 and 1 */
	fReal rand(const Eigen::Matrix<fReal, 2, 1> vecA) const;

	/*map to spherical coordinates*/
	void mapPToSphere(Eigen::Matrix<float, 3, 1>& pos) const;
	void mapVToSphere(Eigen::Matrix<float, 3, 1>& pos, Eigen::Matrix<float, 3, 1>& vel) const;
	/*map to cylindrical coordinates*/
	void mapToCylinder(Eigen::Matrix<float, 3, 1>& pos) const;

	/* Duplicate of quantity's get index */
	inline size_t getIndex(size_t x, size_t y);

public:
	KaminoSolver(size_t nx, size_t ny, fReal radius, fReal gridLength, fReal frameDuration = 1.0 / 30.0);
	~KaminoSolver();

	void stepForward(fReal timeStep);

	void addAttr(std::string name, fReal xOffset = 0.0, fReal yOffset = 0.0);

	void precomputeLaplacian();
	
	KaminoQuantity* getAttributeNamed(std::string name);
	KaminoQuantity* operator[](std::string name);

	void write_data_bgeo(const std::string& s, const int frame);
};