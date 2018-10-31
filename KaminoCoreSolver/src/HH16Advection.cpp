# include "../include/HH16Quantity.h"


void HH16Solver::advectAttrAt(HH16Quantity* attr, size_t gridPhi, size_t gridTheta)
{
    HH16Quantity* uPhi = (*this)["u"];
    HH16Quantity* uTheta = (*this)["v"];

    fReal gPhi = attr->getPhiCoordAtIndex(gridPhi);
    fReal gTheta = attr->getThetaCoordAtIndex(gridTheta);

    fReal guPhi = uPhi->sampleAt(gPhi, gTheta, this->uNorthP, this->uSouthP);
    fReal guTheta = uTheta->sampleAt(gPhi, gTheta, this->uNorthP, this->uSouthP);

    fReal latRadius = this->radius * std::sin(gTheta);
    fReal cofPhi = timeStep / latRadius;
    fReal cofTheta = timeStep / radius;

    fReal deltaPhi = guPhi * cofPhi;
    fReal deltaTheta = guTheta * cofTheta;

    fReal midPhi = gPhi - 0.5 * deltaPhi;
    fReal midTheta = gTheta - 0.5 * deltaTheta;

    fReal muPhi = uPhi->sampleAt(midPhi, midTheta, this->uNorthP, this->uSouthP);
    fReal muTheta = uTheta->sampleAt(midPhi, midTheta, this->uNorthP, this->uSouthP);

    fReal averuPhi = 0.5 * (muPhi + guPhi);
    fReal averuTheta = 0.5 * (muTheta + guTheta);

    fReal newDeltaPhi = averuPhi * cofPhi;
    fReal newDeltaTheta = averuTheta * cofTheta;

    fReal pPhi = gPhi - newDeltaPhi;
    fReal pTheta = gTheta - newDeltaTheta;

    fReal advectedVal = attr->sampleAt(pPhi, pTheta, this->uNorthP, this->uSouthP);
    
    attr->writeValueTo(gridPhi, gridTheta, advectedVal);
}

void HH16Solver::advection()
{
    for (auto quantity : this->attr)
    {
        HH16Quantity* scalarAttr = quantity.second;

        for (size_t gridTheta = 1; gridTheta < scalarAttr->getNTheta() - 1; ++gridTheta)
        {
            for (size_t gridPhi = 0; gridPhi < scalarAttr->getNPhi(); ++gridPhi)
            {
                advectAttrAt(scalarAttr, gridPhi, gridTheta);
            }
        }
    }

	solvePolarVelocitiesAdvection();
	solvePolarScalarsAdvection();

	this->swapAttrBuffers();
}

int signCheck(fReal delta_u, fReal a, fReal b) {
	int sign;
	if (a >= 0 && b >= 0) {
		sign = delta_u >= 0 ? 1 : -1;
	}
	else if (a < 0 && b >= 0) {
		sign = delta_u >= 0 ? -1 : 1;
	}
	else if (a >= 0 && b < 0) {
		sign = delta_u >= 0 ? 1 : -1;
	}
	else {
		sign = delta_u >= 0 ? -1 : 1;
	}
	return sign;
}

// ?? PROBLEM ??
void HH16Solver::solvePolarVelocitiesAdvection()
{
	HH16Quantity* uPhi = (*this)["u"];
	HH16Quantity* uTheta = (*this)["v"];
	fReal u_n, u_down, u_up, u_star, delta_u;
	int sign;

	// HH16 central differencing scheme

	for (size_t gridPhi = 0; gridPhi; ++gridPhi)
	{
		size_t gridShift = (gridPhi + nPhi / 2) % nPhi;

		// North pole
		u_n = uTheta->getValueAt(gridPhi, 0);
		u_down = uTheta->getValueAt(gridPhi, 1);
		u_up = uTheta->getValueAt(gridShift, 1);
		delta_u = u_down - u_n;
		//delta_u = abs(u_up) - abs(u_down);
		//sign = signCheck(delta_u, u_up, u_down);
		//delta_u *= sign;
		//delta_u = u_up - u_down;
		u_star = u_n - timeStep * (u_n / radius) * (delta_u / (1.0 * gridLen));
		this->NPBuffer[gridPhi] = u_star;
		//this->NPBuffer[gridShift] = -u_star;

		// South pole
		u_n = uTheta->getValueAt(gridPhi, nTheta - 1);
		u_down = uTheta->getValueAt(gridPhi, nTheta - 2);
		u_up = uTheta->getValueAt(gridShift, nTheta - 2);
		delta_u = u_n - u_down;
		//delta_u = abs(u_up) - abs(u_down);
		//sign = signCheck(delta_u, u_up, u_down); // South pole has a flipped orientation
		//delta_u *= sign;
		//delta_u = u_up - u_down;
		u_star = u_n - timeStep * (u_n / radius) * (delta_u / (1.0 * gridLen));
		this->SPBuffer[gridPhi] = u_star;
		//this->SPBuffer[gridShift] = -u_star;
	}

	applyPolarBoundaryCondition();
}

// ?? PROBLEM ??
void HH16Solver::solvePolarScalarsAdvection()
{
	// HH16 central differencing scheme

	HH16Quantity* d = (*this)["density"];

	fReal phi_NP, phi_SP;
	size_t phi_NP_idx, phi_SP_idx, gridShift;
	fReal d_n, d_down, d_up, coeff, d_star_NP, d_star_SP;
	fReal delta_d;
	int sign;

	// North pole
	phi_NP = atan2(uNorthP[1], uNorthP[0]); // [-PI, PI]
	phi_NP = phi_NP >= 0 ? phi_NP : M_2PI + phi_NP; // [0, 2PI]
	phi_NP *= invGridLen; // [0, NPhi]
	phi_NP_idx = std::floor(phi_NP);
	phi_NP_idx = phi_NP_idx == nPhi ? 0 : phi_NP_idx; // boundary condition
	gridShift = (phi_NP_idx + nPhi / 2) % nPhi; // opposite polar point

	d_n = d->getValueAt(phi_NP_idx, 0);
	d_down = d->getValueAt(phi_NP_idx, 1);
	d_up = d->getValueAt(gridShift, 1);
	delta_d = abs(d_up) - abs(d_down);
	sign = signCheck(delta_d, d_up, d_down);
	delta_d *= sign;

	coeff = sqrt(uNorthP[0] * uNorthP[0] + uNorthP[1] * uNorthP[1]);

	d_star_NP = d_n - timeStep * (coeff / radius) * (delta_d / (2.0 * gridLen));

	// South pole
	phi_SP = atan2(uSouthP[1], uSouthP[0]);  // [-PI, PI]
	phi_SP = phi_SP >= 0 ? phi_SP : M_2PI + phi_SP;  // [0, 2PI]
	phi_SP *= invGridLen;  // [0, NPhi]
	phi_SP_idx = std::floor(phi_SP);
	phi_SP_idx = phi_SP_idx == nPhi ? 0 : phi_SP_idx; // boundary condition
	gridShift = (phi_SP_idx + nPhi / 2) % nPhi;  // opposite polar point

	d_n = d->getValueAt(phi_SP_idx, nTheta - 1);
	d_down = d->getValueAt(phi_SP_idx, nTheta - 2);
	d_up = d->getValueAt(gridShift, nTheta - 2);
	delta_d = abs(d_up) - abs(d_down);
	sign = -signCheck(delta_d, d_up, d_down); // south pole has a flipped orientation
	delta_d *= sign;
	
	coeff = sqrt(uSouthP[0] * uSouthP[0] + uSouthP[1] * uSouthP[1]);
	
	d_star_SP = d_n - timeStep * (coeff / radius) * (delta_d / (2.0 * gridLen));

	// all pole grid points get same density value
	for (size_t gridPhi = 0; gridPhi < nPhi; ++gridPhi)
	{
		d->writeValueTo(gridPhi, 0, d_star_NP);
		d->writeValueTo(gridPhi, nTheta - 1, d_star_SP);
	}
}

void HH16Solver::resetPoleVelocities()
{
    for (unsigned i = 0; i < 2; ++i)
    {
        uNorthP[i] = 0.0;
        uSouthP[i] = 0.0;
    }
}
