# include "../include/KaminoQuantity.h"

void KaminoSolver::initialize_pressure()
{
	for (size_t i = 0; i < nPhi; ++i) {
		for (size_t j = 0; j < nTheta; ++j) {
			centeredAttr["p"]->setValueAt(i, j, 0.0);
		}
	}
}

void KaminoSolver::initialize_velocity()
{
	KaminoQuantity* u = this->staggeredAttr["u"];
	KaminoQuantity* v = this->staggeredAttr["v"];

	fReal gain = 4096.0 / nPhi;

	for (size_t j = 0; j < u->getNTheta(); ++j)
	{
		for (size_t i = 1; i < u->getNPhi(); ++i)
		{
			fReal ur_x = i * gridLen + gridLen / 2;
			fReal ur_y = (j + 1) * gridLen;
			fReal lr_x = i * gridLen + gridLen / 2;
			fReal lr_y = j * gridLen;
			fReal ul_x = i * gridLen - gridLen / 2;
			fReal ul_y = (j + 1) * gridLen;
			fReal ll_x = i * gridLen - gridLen / 2;
			fReal ll_y = j * gridLen;
			fReal noise_ur = FBM(ur_x, ur_y);
			fReal noise_lr = FBM(lr_x, lr_y);
			fReal noise_ul = FBM(ul_x, ul_y);
			fReal noise_ll = FBM(ll_x, ll_y);
			fReal noiseDy_l = (noise_ur - noise_lr) / (radius * gridLen);
			fReal noiseDy_r = (noise_ul - noise_ll) / (radius * gridLen);
			fReal avgNoise = (noiseDy_l + noiseDy_r) / 2.0;
			u->setValueAt(i, j, avgNoise * gain);
		}
	}
	// phi = 0 seam
	for (size_t j = 0; j < u->getNTheta(); ++j)
	{
		fReal ur_x = gridLen / 2;
		fReal ur_y = (j + 1) * gridLen;
		fReal lr_x = gridLen / 2;
		fReal lr_y = j * gridLen;
		fReal ul_x = 2 * M_PI - gridLen / 2;
		fReal ul_y = (j + 1) * gridLen;
		fReal ll_x = 2 * M_PI - gridLen / 2;
		fReal ll_y = j * gridLen;
		fReal noise_ur = FBM(ur_x, ur_y);
		fReal noise_lr = FBM(lr_x, lr_y);
		fReal noise_ul = FBM(ul_x, ul_y);
		fReal noise_ll = FBM(ll_x, ll_y);
		fReal noiseDy_l = (noise_ur - noise_lr) / (radius * gridLen);
		fReal noiseDy_r = (noise_ul - noise_ll) / (radius * gridLen);
		fReal avgNoise = (noiseDy_l + noiseDy_r) / 2.0;
		u->setValueAt(0, j, avgNoise * gain);
	}

	// u_theta at poles is set to zero
	for (size_t i = 0; i < v->getNPhi(); ++i)
	{
		v->setValueAt(i, 0, 0);
		v->setValueAt(i, v->getNTheta() - 1, 0);
	}

	// set u_theta initial values using FBM curl noise
	for (size_t j = 1; j < v->getNTheta() - 1; ++j)
	{
		for (size_t i = 0; i < v->getNPhi(); ++i)
		{
			fReal ur_x = (i + 1) * gridLen;
			fReal ur_y = j * gridLen + gridLen / 2;
			fReal lr_x = (i + 1) * gridLen;
			fReal lr_y = j * gridLen - gridLen / 2;
			fReal ul_x = i * gridLen;
			fReal ul_y = j * gridLen + gridLen / 2;
			fReal ll_x = i * gridLen;
			fReal ll_y = j * gridLen + gridLen / 2;
			fReal noise_ur = FBM(ur_x, ur_y);
			fReal noise_lr = FBM(lr_x, lr_y);
			fReal noise_ul = FBM(ul_x, ul_y);
			fReal noise_ll = FBM(ll_x, ll_y);
			fReal noiseDy_u = -1 * (noise_ur - noise_ul) / (radius * gridLen);
			fReal noiseDy_d = -1 * (noise_lr - noise_ll) / (radius * gridLen);
			fReal avgNoise = (noiseDy_u + noiseDy_d) / 2.0;
			v->setValueAt(i, j, avgNoise * gain);
		}
	}

	// Solve the polar velocities first.
	solvePolarVelocities();

	// Set up the skewed values.
	
	size_t skewedPhiInd = 0;
	size_t skewedThetaInd = 0;
	for (int j = 0; j < u->getNTheta(); ++j)
	{
		for (int i = 0; i < u->getNPhi(); ++i)
		{
			u->convert2SlewedCoord(i, j, skewedPhiInd, skewedThetaInd);
			u->writeValueTo(i, j, u->getValueAt(skewedPhiInd, skewedThetaInd));
		}
	}
	for (int j = 1; j < v->getNTheta() - 1; ++j)
	{
		for (int i = 0; i < v->getNPhi(); ++i)
		{
			v->convert2SlewedCoord(i, j, skewedPhiInd, skewedThetaInd);
			v->writeValueTo(i, j, v->getValueAt(skewedPhiInd, skewedThetaInd));
		}
	}

	// Copy back to this buffer.
	for (size_t j = 0; j < u->getNTheta(); ++j)
	{
		for (size_t i = 0; i < u->getNPhi(); ++i)
		{
			u->setValueAt(i, j, u->getNextValueAt(i, j));
		}
	}

	for (size_t j = 1; j < v->getNTheta() - 1; ++j)
	{
		for (size_t i = 0; i < v->getNPhi(); ++i)
		{
			v->setValueAt(i, j, v->getNextValueAt(i, j));
		}
	}

	solvePolarVelocities();
}

enum componentsSPHERE { radiusComp, phiComp, thetaComp };
enum componentsXYZ { xComp, yComp, zComp };

Eigen::Vector3d spherical2Cartesian(const Eigen::Vector3d& input)
{
	Eigen::Vector3d ret = Eigen::Vector3d::Zero();
	ret[xComp] = input[radiusComp] * std::sin(input[thetaComp]) * std::cos(input[phiComp]);
	ret[yComp] = input[radiusComp] * std::sin(input[thetaComp]) * std::sin(input[phiComp]);
	ret[zComp] = input[radiusComp] * std::cos(input[thetaComp]);

	return ret;
}

Eigen::Vector3d sphericalCrossProd(const Eigen::Vector3d& omega, const Eigen::Vector3d& r)
{
	//Omega and r are in spherical coordinates. Component order: see the enum above
	Eigen::Vector3d omegaXyz = spherical2Cartesian(omega);
	Eigen::Vector3d rXyz = spherical2Cartesian(r);
	Eigen::Vector3d vel = omegaXyz.cross(rXyz);

	fReal phi = r[phiComp];
	fReal theta = r[thetaComp];
	fReal vx = vel[xComp];
	fReal vy = vel[yComp];
	fReal vz = vel[zComp];

	fReal vPhi = -vx * std::sin(phi) + vy * std::cos(phi);
	fReal vProj = vx * std::cos(phi) + vy * std::sin(phi);

	fReal vTheta = vProj * std::cos(theta) - vz * std::sin(theta);
	fReal vR = vProj * std::sin(theta) + vz * std::cos(theta);

	Eigen::Vector3d ret = Eigen::Vector3d::Zero();
	ret[radiusComp] = vR;
	ret[phiComp] = vPhi;
	ret[thetaComp] = vTheta;

	return ret;
}

void KaminoSolver::initializeVelocityFromOmega(Eigen::Vector3d omega)
{
	KaminoQuantity* u = this->staggeredAttr["u"];
	KaminoQuantity* v = this->staggeredAttr["v"];

	for (size_t beltJ = 0; beltJ < this->nTheta; ++beltJ)
	{
		fReal beltTheta = (static_cast<fReal>(beltJ) + 0.5) * gridLen;
		for (size_t gridI = 0; gridI < this->nPhi; ++gridI)
		{
			fReal beltPhi = (static_cast<fReal>(gridI)) * gridLen;

			fReal phiLeft = gridI == 0 ? M_2PI - 0.5 * gridLen : beltPhi - 0.5 * gridLen;
			Eigen::Vector3d r(radius, phiLeft, beltTheta);
			Eigen::Vector3d vel = sphericalCrossProd(omega, r);
			u->setValueAt(gridI, beltJ, vel[phiComp]);

			fReal phiRight = beltPhi + 0.5 * gridLen;
			r = Eigen::Vector3d(radius, phiRight, beltTheta);
			vel = sphericalCrossProd(omega, r);
			u->setValueAt((gridI + 1) % nPhi, beltJ, vel[phiComp]);

			fReal thetaLower = beltTheta - 0.5 * gridLen;
			r = Eigen::Vector3d(radius, beltPhi, thetaLower);
			vel = sphericalCrossProd(omega, r);
			v->setValueAt(gridI, beltJ, vel[thetaComp]);

			fReal thetaHigher = beltTheta + 0.5 * gridLen;
			r = Eigen::Vector3d(radius, beltPhi, thetaHigher);
			vel = sphericalCrossProd(omega, r);
			v->setValueAt(gridI, beltJ + 1, vel[thetaComp]);
		}
	}

	// Heat up the next buffer.
	for (size_t j = 0; j < u->getNTheta(); ++j)
	{
		for (size_t i = 0; i < u->getNPhi(); ++i)
		{
			u->writeValueTo(i, j, u->getValueAt(i, j));
		}
	}
	for (size_t j = 0; j < v->getNTheta(); ++j)
	{
		for (size_t i = 0; i < v->getNPhi(); ++i)
		{
			v->writeValueTo(i, j, v->getValueAt(i, j));
		}
	}

	solvePolarVelocities();
	u->swapBuffer();
	v->swapBuffer();
}

fReal KaminoSolver::fPhi(const fReal x)
{
	fReal arg = x;
	return sin(arg) + 0.1 * B * std::cos(std::rand());
}

fReal KaminoSolver::gTheta(const fReal y)
{
	fReal arg = y;
	return cos(arg) + 0.1 * C * std::sin(std::rand());
}

fReal KaminoSolver::lPhi(const fReal x)
{
	fReal arg = x;
	return cos(arg) + 0.1 * D * std::sin(std::rand());
}

fReal KaminoSolver::mTheta(const fReal y)
{
	fReal arg = y;
	return sin(arg) + 0.1 * E * std::cos(std::rand());
}

fReal KaminoSolver::FBM(const fReal x, const fReal y) {
	fReal total = 0.0f;
	fReal resolutionX = 0.15;
	fReal resolutionY = 0.5;
	fReal persistance = 0.5;
	int octaves = 4;

	for (int i = 0; i < octaves; i++) {
		fReal freq = std::pow(2.0f, i);
		fReal amp = std::pow(persistance, i);
		total += amp * interpNoise2D(x * freq / resolutionX, y * freq / resolutionY);
	}
	fReal a = 1 - persistance;  // normalization

	return a * total / 2.0f;  // normalized, pseudorandom number between -1 and 1
}

fReal KaminoSolver::interpNoise2D(const fReal x, const fReal y) const {
	fReal intX = std::floor(x);
	fReal fractX = x - intX;
	fReal intY = std::floor(y);
	fReal fractY = y - intY;

	fReal v1 = rand(Eigen::Matrix<fReal, 2, 1>(intX, intY));
	fReal v2 = rand(Eigen::Matrix<fReal, 2, 1>(intX + 1, intY));
	fReal v3 = rand(Eigen::Matrix<fReal, 2, 1>(intX, intY + 1));
	fReal v4 = rand(Eigen::Matrix<fReal, 2, 1>(intX + 1, intY + 1));

	// interpolate for smooth transitions
	fReal i1 = KaminoLerp(v1, v2, fractX);
	fReal i2 = KaminoLerp(v3, v4, fractX);
	return KaminoLerp(i1, i2, fractY);
}

fReal KaminoSolver::rand(const Eigen::Matrix<fReal, 2, 1> vecA) const {
	// return pseudorandom number between -1 and 1
	Eigen::Matrix<fReal, 2, 1> vecB = Eigen::Matrix<fReal, 2, 1>(12.9898, 4.1414);
	fReal val = sin(vecA.dot(vecB) * 43758.5453);
	return val - std::floor(val);
}

void KaminoSolver::initialize_density()
{
	for (size_t i = 0; i < nPhi; ++i)
	{
		for (size_t j = 0; j < nTheta; ++j)
		{
			centeredAttr["density"]->setValueAt(i, j, 0.0);
		}
	}
}

void KaminoSolver::initialize_boundary()
{
	for (size_t gridX = 0; gridX != this->nPhi / 2; ++gridX)
	{
		//this->gridTypes[getIndex(gridX, nTheta / 2)] = SOLIDGRID;
		//this->gridTypes[getIndex(gridX, this->nTheta - 1)] = SOLIDGRID;
	}
	/*for (size_t gridY = 0; gridY != this->nTheta; ++gridY)
	{
		this->gridTypes[getIndex(0, gridY)] = SOLIDGRID;
		this->gridTypes[getIndex(nPhi / 2, gridY)] = SOLIDGRID;
	}*/
}
