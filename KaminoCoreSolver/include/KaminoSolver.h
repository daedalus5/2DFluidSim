# pragma once

# include <string>
# include <map>
# include <iostream>
# include <vector>
# include "Partio.h"
# include <Eigen/Core>
# include <Eigen/Dense>
# include <cmath>

# define DEBUGBUILD

// The solution to switch between double and float
typedef double fReal;


// Handy Lerp.
template <class Type>
Type KaminoLerp(const Type &fromEndPoint, const Type &toEndPoint, double factor)
{
	return (1.0 - factor) * fromEndPoint + factor * toEndPoint;
}


// The attribute base class.
class KaminoAttribute
{
private:
	
	std::string attrName;

protected:
	/* Grid dimensions */
	/* These are only assigned by the Grid and not mutable */
	size_t nx;
	size_t ny;

	/* Grid size */
	fReal gridLen;

	/* Double buffer */
	fReal* thisStep;
	fReal* nextStep;

	/* Swap the buffer */
	void swapBuffer();

public:
	/* Constructor */
	KaminoAttribute(std::string attributeName, size_t nx, size_t ny, fReal gridLen);
	/* Destructor */
	virtual ~KaminoAttribute();

	/* Getter */
	virtual fReal getValueAt(size_t x, size_t y);
	/* Setter */
	virtual void setValueAt(size_t x, size_t y, fReal val);
	/* Access */
	virtual fReal& accessValueAt(size_t x, size_t y) = 0;
	/* Lerped Sampler using world coordinates */
	virtual fReal sampleAt(fReal x, fReal y);
	/* Lerped Sampler taking in grid coordinates (treat gridLen as 1.0) */
	//virtual fReal sampleAtGC(fReal x, fReal y) = 0;
};


// Those that are centered.
/*	
	Index conventions: 
	All indices start at zero, with the centered attribute of the first grid as attribute(0, 0)
	The corresponding faced attribute to the left is u0, and to the right is u1
	That goes the same with v: to the lower y is v0, and to the higher is v1
*/
class KaminoCenteredAttr : public KaminoAttribute
{
public:
	KaminoCenteredAttr(std::string attributeName, size_t nx, size_t ny, fReal gridLen);
	virtual ~KaminoCenteredAttr();

	/* Access */
	fReal& accessValueAt(size_t x, size_t y) override;
	/* Lerped Sampler */
	//fReal sampleAtGC(fReal x, fReal y) override;
};


// The U velocity.
/*
	Note that we have to allocate for (nx + 1) by ny fReals
*/
class KaminoUAttr : public KaminoAttribute
{
public:
	KaminoUAttr(std::string attributeName, size_t nx, size_t ny, fReal gridLen);
	virtual ~KaminoUAttr();

	/* Access */
	fReal& accessValueAt(size_t x, size_t y) override;
	/* Lerped Sampler */
	//fReal sampleAtGC(fReal x, fReal y) override;
};


// The V velocity
/*
	Note that we have to allocate for nx by (ny + 1) fReals
*/
class KaminoVAttr : public KaminoAttribute
{
public:
	KaminoVAttr(std::string attributeName, size_t nx, size_t ny, fReal gridLen);
	virtual ~KaminoVAttr();

	/* Access */
	fReal& accessValueAt(size_t x, size_t y) override;
	/* Lerped Sampler */
	//fReal sampleAtGC(fReal x, fReal y) override;
};

// The solver class.
class KaminoGrid
{
private:
	/* Grid dimensions */
	size_t nx;
	size_t ny;
	/* Grid size */
	fReal gridLen;

	/* So that it remembers all these attributes within */
	std::map<std::string, KaminoAttribute*> attributeTable;

	/* Something about time steps */
	fReal frameDuration;
	fReal timeStep;
	fReal timeElapsed;

	// TODO
	// void advection();
	// void projection();
	// void bodyForce();

	/* distribute initial velocity values at grid points */
    void initialize_velocity();
    /* FBM noise function for velocity distribution */
    fReal FBM(const fReal x, const fReal y, const fReal persistance, const int octaves) const;
    /* 2D noise interpolation function for smooth FBM noise */
    fReal interpNoise2D(const fReal x, const fReal y) const;
    /* returns a pseudorandom number between -1 and 1 */
    fReal rand(const Eigen::Matrix<fReal, 2, 1> vecA) const;

public:
	KaminoGrid(size_t nx, size_t ny, fReal gridLength, fReal frameDuration = 1.0 / 30.0);
	~KaminoGrid();

	void stepForward(fReal timeStep);

	void addCenteredAttr(std::string name);
	void addUAttr(std::string name);
	void addVAttr(std::string name);
	KaminoAttribute* getAttributeNamed(std::string name);
	KaminoAttribute* operator[](std::string name);

	void write_data_bgeo(const std::string& s, const int frame);
};