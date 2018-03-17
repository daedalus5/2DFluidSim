# include "../include/KaminoSolver.h"

KaminoGrid::KaminoGrid(size_t nx, size_t ny, fReal gridLength, fReal frameDuration) :
	nx(nx), ny(ny), gridLen(gridLength), frameDuration(frameDuration),
	timeStep(0.0), timeElapsed(0.0)
{
	addUAttr("u");		// u velocity
	addVAttr("v");		// v velocity
	addCenteredAttr("rho");	// rho density
}

KaminoGrid::~KaminoGrid()
{
	for (auto& attr : this->attributeTable)
	{
		delete attr.second;
	}
}

void KaminoGrid::stepForward(fReal timeStep)
{
	this->timeStep = timeStep;
	advection();
	projection();
	bodyForce();
}

void KaminoGrid::addCenteredAttr(std::string name)
{
	KaminoAttribute* ptr = new KaminoCenteredAttr(name, this->nx, this->ny, this->gridLen);
	this->attributeTable.emplace(std::pair<std::string, KaminoAttribute*>(name, ptr));
}

void KaminoGrid::addUAttr(std::string name)
{
	KaminoAttribute* ptr = new KaminoUAttr(name, this->nx, this->ny, this->gridLen);
	this->attributeTable.emplace(std::pair<std::string, KaminoAttribute*>(name, ptr));
}

void KaminoGrid::addVAttr(std::string name)
{
	KaminoAttribute* ptr = new KaminoVAttr(name, this->nx, this->ny, this->gridLen);
	this->attributeTable.emplace(std::pair<std::string, KaminoAttribute*>(name, ptr));
}

KaminoAttribute* KaminoGrid::getAttributeNamed(std::string name)
{
	return (*this).getAttributeNamed(name);
}

KaminoAttribute* KaminoGrid::operator[](std::string name)
{
	return this->attributeTable.at(name);
}

////////////// NEW ////////////////////////
fReal KaminoGrid::FBM(const fReal x, const fReal y, const fReal persistance, const int octaves) const{
    fReal total = 0.0f;
    fReal resolution = 1.f;

    for(int i = 0; i < octaves; i++){
        fReal freq = std::pow(2.0f, i);
        fReal amp = std::pow(persistance, i);
        total += amp * interpNoise2D(x * freq / resolution, y * freq / resolution);
    }
    fReal a = 1 - persistance;  // normalization

    return a * (1 + total) / 2.0f;  // normalized, pseudorandom number between 0 and 1

}

fReal KaminoGrid::interpNoise2D(const fReal x, const fReal y) const{
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

fReal KaminoGrid::rand(const Eigen::Matrix<fReal, 2, 1> vecA) const{
    // return pseudorandom number between -1 and 1
    Eigen::Matrix<fReal, 2, 1> vecB = Eigen::Matrix<fReal, 2, 1>(12.9898, 4.1414);
	fReal val = sin(vecA.dot(vecB) * 43758.5453);
	return val - std::floor(val);
}


///////////// FIX //////////////////////

void KaminoGrid::write_data_bgeo(const std::string& s, const int frame)
{
# ifndef _MSC_VER
    std::string file = s + std::to_string(frame) + ".bgeo";

    Eigen::Matrix<float, 2, 1> gridPos[nx][ny];
    Eigen::Matrix<float, 2, 1> gridVel[nx][ny];

    float x = gridLen / 2.0;
    float y = gridLen / 2.0;
    for(size_t i = 0; i < nx; ++i){
    	for(size_t j = 0; j < ny; ++j){
    		gridPos[i][j] = Eigen::Matrix<float, 2, 1>(x, y);
    		gridVel[i][j] = Eigen::Matrix<float, 2, 1>(0.0, 0.0);
    		x += gridLen;
    	}
    	y += gridLen;
    }
    // TODO: interpolate velocities to grid centers and combine into vec2

    Partio::ParticlesDataMutable* parts = Partio::create();
    Partio::ParticleAttribute pH, vH;
    pH = parts->addAttribute("p", Partio::VECTOR, 2);
    vH = parts->addAttribute("v", Partio::VECTOR, 2);

    for(size_t i = 0; i < ny; ++i){
        for(size_t j = 0; j < nx; ++j){
            int idx = parts->addParticle();
            float* p = parts->dataWrite<float>(pH, idx);
            float* v = parts->dataWrite<float>(vH, idx);
            for (int k = 0; k < 2; ++k){
                p[k] = gridPos[i][j][k];
                v[k] = gridVel[i][j][k];
            }
        }
    }
    Partio::write(file.c_str(), *parts);
    parts->release();
# endif
}

void KaminoGrid::build_particle_grid()
{
	// fReal x = gridLen / 2.0;
	// fReal y = gridLen / 2.0;
 //    for(size_t i = 0; i < nx; ++i){
 //    	for(size_t j = 0; j < ny; ++j){
 //    		//
 //    	}
 //    }
}

void KaminoGrid::distribute_velocity()
{
    // for(int i = 0; i < nx; ++i){
    //     for(int j = 0; j < ny; ++j){
    //         velocities[i][j] = Eigen::Matrix<float, 2, 1>(0.0, 0.0);
    //     }
    // }
}