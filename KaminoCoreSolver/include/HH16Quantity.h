# pragma once

# include "GlobalIncludes.h"
# include <string>
# include <map>
# include <iostream>
# include <vector>
//# ifndef _MSC_VER
# include "Partio.h"
//# endif
# include <Eigen/Core>
# include <Eigen/Dense>
# include <cmath>

# define OMParallelize

# ifdef OMParallelize
# include <omp.h>
# define TOTALThreads 16
# endif

# define DEBUGBUILD

//const fReal density = 1000.0;

// Handy Lerp.
template <class Type>
Type Lerp(const Type &fromEndPoint, const Type &toEndPoint, double factor)
{
    return (1.0 - factor) * fromEndPoint + factor * toEndPoint;
}

//enum Coord { x, y };

// The attribute base class.
class HH16Quantity
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

    /* Double buffer */
    fReal* thisStep;
    fReal* nextStep;

    /* Get index */
    size_t getIndex(size_t x, size_t y);

public:
    /* Constructor */
    HH16Quantity(std::string attributeName, size_t nx, size_t ny, fReal gridLen);
    /* Destructor */
    ~HH16Quantity();

    /* Swap the buffer */
    void swapBuffer();
    /* Get its name */
    std::string getName();
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
    fReal sampleAt(fReal x, fReal y, fReal uNorth[2], fReal uSouth[2]);
    /* Given the index, show its origin in world coordinates*/
    fReal getPhiCoordAtIndex(size_t phi);
    fReal getThetaCoordAtIndex(size_t theta);
    /* And given world coordinates, show its index backwards... */
    size_t getPhiIndexAtCoord(fReal phi);
    size_t getThetaIndexAtCoord(fReal theta);
};

bool validatePhiTheta(fReal & phi, fReal & theta);

// The solver class.
class HH16Solver
{
private:
    /* Grid dimensions */
    size_t nPhi;
    size_t nTheta;
    /* Radius of sphere */
    fReal radius;
    /* Grid size */
    fReal gridLen;
    /* Inverted grid size*/
    fReal invGridLen;

    /* So that it remembers all these attributes within */
    std::map<std::string, HH16Quantity*> attr;

    /* Something about time steps */
    fReal frameDuration;
    fReal timeStep;
    fReal timeElapsed;

    float advectionTime;
    float geometricTime;
    float projectionTime;

    void resetPoleVelocities();
    void averageVelocities();
    void solvePolarVelocities();

    // We only have to treat uTheta differently
    void advectAttrAt(HH16Quantity* attr, size_t gridPhi, size_t gridTheta);

    void advectionScalar();
    void advectionSpeed();

    void geometric();
    void projection();
    void bodyForce();

    // Run Intel MKL's Spherical Poisson Solver
    void MKLSphericalPoisson();

    // Swap all these buffers of the attributes.
    void swapAttrBuffers();

    /* distribute initial velocity values at grid points */
    void initialize_velocity();
    /* initialize pressure attribute */
    void initialize_pressure();
    /* initialize density distribution */
    void initialize_density();

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
    
public:
    // Velocities at poles in xyz cartesian coordinates
    fReal uNorthP[2];
    fReal uSouthP[2];

    HH16Solver(size_t nx, size_t ny, fReal radius, fReal gridLength, fReal frameDuration);
    ~HH16Solver();

    void stepForward(fReal timeStep);

    void addAttr(std::string name);

    HH16Quantity* getAttributeNamed(std::string name);
    HH16Quantity* operator[](std::string name);

    void write_data_bgeo(const std::string& s, const int frame);

    /* Duplicate of quantity's get index */
    size_t getIndex(size_t x, size_t y);
};