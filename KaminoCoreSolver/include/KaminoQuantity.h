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
# include <unsupported/Eigen/FFT>

# define M_PI           3.14159265358979323846  /* pi */
# define M_2PI			6.28318530717958647692  /* 2pi */
# define M_hPI			1.57079632679489661923  /* pi / 2*/

# define DEBUGBUILD

// The solution to switch between double and float
typedef double fReal;

const fReal density = 1000.0;
const fReal uSolid = 0.0;
const fReal vSolid = 0.0;

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
	fReal getPhiCoordAtIndex(size_t phi);
	fReal getThetaCoordAtIndex(size_t theta);
	/* And given world coordinates, show its index backwards... */
	size_t getPhiIndexAtCoord(fReal phi);
	size_t getThetaIndexAtCoord(fReal theta);

	fReal getPhiOffset();
	fReal getThetaOffset();
};

void validatePhiTheta(fReal & phi, fReal & theta);

struct tracer
{
	fReal phi;
	fReal theta;
	fReal radius;
	tracer(fReal phi, fReal theta, fReal radius) : phi(phi), theta(theta), radius(radius)
	{}
	void tracerStepForward(fReal uPhi, fReal uTheta, fReal timeStep)
	{
		this->phi += timeStep * uPhi / (radius);
		this->theta += timeStep * uTheta / (radius * std::sin(theta));
		validatePhiTheta(phi, theta);
	}
	void getCartesianXYZ(fReal& x, fReal& y, fReal& z)
	{
		x = radius * std::sin(theta) * std::cos(phi);
		z = radius * std::sin(theta) * std::sin(phi);
		y = radius * std::cos(theta);
	}
};

// The solver class.
class KaminoSolver
{
private:
	// Buffer for the capital U.
	fReal* fourierU;
	// Buffer for the divergence, before the transform.
	fReal* beffourierF;
	// Buffer for the divergence, F n theta.
	fReal* fourieredF;
	// Diagonal elements a (lower);
	fReal* a;
	// Diagonal elements b (major diagonal);
	fReal* b;
	// Diagonal elements c (upper);
	fReal* c;
	// Divergence fourier coefficients
	fReal* d;

	/* Grid types */
	gridType* gridTypes;
	/* Grid dimensions */
	size_t nPhi;
	size_t nTheta;
	/* Radius of sphere */
	fReal radius;
	/* Grid size */
	fReal gridLen;
	/* Inverted grid size*/
	fReal invGridLen;
	/* Laplacian Matrix */
	//Eigen::SparseMatrix<fReal> Laplacian;

	/* So that it remembers all these attributes within */
	std::map<std::string, KaminoQuantity*> centeredAttr;
	std::map<std::string, KaminoQuantity*> staggeredAttr;

	/* Something about time steps */
	fReal frameDuration;
	fReal timeStep;
	fReal timeElapsed;

	// Velocities at poles in xyz cartesian coordinates
	//fReal uThetaNorthP[2];
	fReal uPhiNorthP[2];
	//fReal uThetaSouthP[2];
	fReal uPhiSouthP[2];

	tracer trc;

	void resetPoleVelocities();
	void averageVelocities();

	// Is it solid? or fluid? or even air?
	gridType getGridTypeAt(size_t x, size_t y);

	// We only have to treat uTheta differently
	void advectAttrAt(KaminoQuantity* attr, size_t gridPhi, size_t gridTheta);

	void advectionScalar();
	void advectionSpeed();

	void geometric();
	void projection();
	void bodyForce();
	void updateTracer();

	void fillDivergence();
	void transformDivergence();
	void invTransformPressure();

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

	/* Tri-diagonal matrix solver */
	void TDMSolve(fReal* a, fReal* b, fReal* c, fReal* d);
	/* Load diagonal element arrays */
	void loadABC(size_t n);

public:
	KaminoSolver(size_t nx, size_t ny, fReal radius, fReal gridLength, fReal frameDuration = 1.0 / 30.0);
	~KaminoSolver();

	void stepForward(fReal timeStep);

	void addCenteredAttr(std::string name, fReal xOffset = 0.5, fReal yOffset = 0.5);
	void addStaggeredAttr(std::string name, fReal xOffset, fReal yOffset);

	//void precomputeLaplacian();
	
	KaminoQuantity* getAttributeNamed(std::string name);
	KaminoQuantity* operator[](std::string name);

	void write_data_bgeo(const std::string& s, const int frame);
	void write_data_tracer(const std::string& s, const int frame);
};