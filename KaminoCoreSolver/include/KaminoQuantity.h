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
	KaminoQuantity(std::string attributeName, size_t nx, size_t ny, fReal gridLen);
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
