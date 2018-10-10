# include "../include/HH16Quantity.h"

void HH16Solver::initialize_velocity()
{
    HH16Quantity* u = this->attr["u"];
    HH16Quantity* v = this->attr["v"];

    // set u_phi initial values using FBM curl noise
    fReal gain = 4096.0 / nPhi;

	for (size_t i = 0; i < u->getNPhi(); ++i) {
		for (size_t j = 0; j < u->getNTheta(); ++j) {
			//fReal val = FBM(cos(i * gridLen), sin(j * gridLen));
			fReal val = sin(i * gridLen) * sin(j * gridLen);
			u->setValueAt(i, j, val);
		}
	}

	for (size_t i = 0; i < v->getNPhi(); ++i) {
		for (size_t j = 0; j < v->getNTheta(); ++j) {
			//fReal val = FBM(cos(i * gridLen), sin(j * gridLen));
			fReal val = sin(i * gridLen) * sin(j * gridLen);
			v->setValueAt(i, j, val);
		}
	}

	/*
    for (size_t j = 0; j < u->getNTheta() - 1; ++j)
    {
        for (size_t i = 1; i < u->getNPhi(); ++i)
        {
            fReal ur_x = i * gridLen;
            fReal ur_y = (j + 1) * gridLen;
            fReal lr_x = i * gridLen;
            fReal lr_y = j * gridLen;
            fReal ul_x = (i - 1) * gridLen;
            fReal ul_y = (j + 1) * gridLen;
            fReal ll_x = (i - 1) * gridLen;
            fReal ll_y = j * gridLen;
            fReal noise_ur = FBM(sin(ur_x), ur_y);
            fReal noise_lr = FBM(sin(lr_x), lr_y);
            fReal noise_ul = FBM(sin(ul_x), ul_y);
            fReal noise_ll = FBM(sin(ll_x), ll_y);
            fReal noiseDy_l = (noise_ur - noise_lr) / (radius * gridLen);
            fReal noiseDy_r = (noise_ul - noise_ll) / (radius * gridLen);
            fReal avgNoise = (noiseDy_l + noiseDy_r) / 2.0;
            u->setValueAt(i, j, avgNoise * gain);
        }
    }
    // 0 index for u_phi
    for (size_t j = 0; j < u->getNTheta() - 1; ++j){
        fReal ur_x = 0;
        fReal ur_y = (j + 1) * gridLen;
        fReal lr_x = 0;
        fReal lr_y = j * gridLen;
        fReal ul_x = 2 * M_PI - gridLen;
        fReal ul_y = (j + 1) * gridLen;
        fReal ll_x = 2 * M_PI - gridLen;
        fReal ll_y = j * gridLen;
        fReal noise_ur = FBM(sin(ur_x), ur_y);
        fReal noise_lr = FBM(sin(lr_x), lr_y);
        fReal noise_ul = FBM(sin(ul_x), ul_y);
        fReal noise_ll = FBM(sin(ll_x), ll_y);
        fReal noiseDy_l = (noise_ur - noise_lr) / (radius * gridLen);
        fReal noiseDy_r = (noise_ul - noise_ll) / (radius * gridLen);
        fReal avgNoise = (noiseDy_l + noiseDy_r) / 2.0;
        u->setValueAt(0, j, avgNoise * gain);
    }

    // u_theta at poles is set to zero
	// NOTE : could be problematic...
    for(size_t i = 0; i < v->getNPhi(); ++i)
    {
        v->setValueAt(i, 0, 0);
        v->setValueAt(i, v->getNTheta() - 1, 0);
    }

    // set u_theta initial values using FBM curl noise
    for(size_t j = 1; j < v->getNTheta() - 1; ++j)
    {
        for(size_t i = 0 ; i < v->getNPhi(); ++i)
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

	*/

	// should run a projection step here...

	// Heat up the next buffers.
	// Polar values are written in applyPolarBoundaryCondition()

	for (size_t i = 0; i < u->getNPhi(); ++i) {
		for (size_t j = 1; j < u->getNTheta() - 1; ++j) {
			u->writeValueTo(i, j, u->getValueAt(i, j));
			v->writeValueTo(i, j, v->getValueAt(i, j));
		}
	}

    // Heat up the polar buffers.

	for (size_t gridPhi = 0; gridPhi < nPhi; ++gridPhi) {
		this->NPBuffer[gridPhi] = u->getValueAt(gridPhi, 0);
		this->SPBuffer[gridPhi] = u->getValueAt(gridPhi, u->getNTheta() - 1);
	}

	applyPolarBoundaryCondition();

    u->swapBuffer();
    v->swapBuffer();
}

fReal HH16Solver::FBM(const fReal x, const fReal y) {
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

fReal HH16Solver::interpNoise2D(const fReal x, const fReal y) const {
    fReal intX = std::floor(x);
    fReal fractX = x - intX;
    fReal intY = std::floor(y);
    fReal fractY = y - intY;

    fReal v1 = rand(Eigen::Matrix<fReal, 2, 1>(intX, intY));
    fReal v2 = rand(Eigen::Matrix<fReal, 2, 1>(intX + 1, intY));
    fReal v3 = rand(Eigen::Matrix<fReal, 2, 1>(intX, intY + 1));
    fReal v4 = rand(Eigen::Matrix<fReal, 2, 1>(intX + 1, intY + 1));

    // interpolate for smooth transitions
    fReal i1 = Lerp(v1, v2, fractX);
    fReal i2 = Lerp(v3, v4, fractX);
    return Lerp(i1, i2, fractY);
}

fReal HH16Solver::rand(const Eigen::Matrix<fReal, 2, 1> vecA) const {
    // return pseudorandom number between -1 and 1
    Eigen::Matrix<fReal, 2, 1> vecB = Eigen::Matrix<fReal, 2, 1>(12.9898, 4.1414);
    fReal val = sin(vecA.dot(vecB) * 43758.5453);
    return val - std::floor(val);
}

void HH16Solver::initialize_density()
{
    for(size_t i = 0; i < nPhi; ++i)
    {
        for(size_t j = 0; j < nTheta / 2; ++j)
        {
            attr["density"]->setValueAt(i, j, 1.0);
        }
		for (size_t j = nTheta / 2; j < nTheta; ++j)
		{
			attr["density"]->setValueAt(i, j, 0.0);
		}
    }
}

void HH16Solver::initialize_pressure()
{
	for (size_t i = 0; i < nPhi; ++i) {
		for (size_t j = 0; j < nTheta; ++j) {
			attr["pressure"]->setValueAt(i, j, 0.0);
		}
	}
}