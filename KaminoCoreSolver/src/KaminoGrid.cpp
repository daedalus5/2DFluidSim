# include "../include/KaminoSolver.h"

KaminoGrid::KaminoGrid(size_t nx, size_t ny, fReal gridLength, fReal frameDuration) :
	nx(nx), ny(ny), gridLen(gridLength), frameDuration(frameDuration),
	timeStep(0.0), timeElapsed(0.0)
{
	addUAttr("u");		// u velocity
	addVAttr("v");		// v velocity
	addCenteredAttr("rho");	// rho density

	initialize_velocity();
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
	// TODO
	// advection();
	// bodyForce();
	// projection();
}

void KaminoGrid::projection()
{
	fReal density = 1.0;	// rest fluid density
	fReal scale = timeStep / density;
	fReal invGridLen = 1 / gridLen;

	// construct the matrix A
	// present construction assumes 2D fluid in every cell and toroidal BCs
	Eigen::Matrix<fReal, nx*ny, nx*ny> A;
	A.setZero();

	// construct A row-by-row
	size_t k = 1;
	Eigen::Matrix<fReal, nx*ny, 1> row;
	for(size_t i = 0; i < nx; ++i){
		for(size_t j = 0; j < ny; ++j){
			row.setZero();
			row(j*nx + i) = 4;
			i + 1 > nx ? row(j*nx) = -1 : row(j*nx + i + 1) = -1;
			i - 1 < 0 ? row(j*nx + nx - 1) = -1 : row(j*nx + i - 1) = -1;
			j + 1 > ny ? row(i) = -1 : row((j + 1)*nx + i) = -1;
			j - 1 < 0 ? row((ny - 1)*nx + i) = -1 : row((j - 1)*nx + i) = -1;
			A.row(k) = row;
			k++;
		}
	}
	A *= scale;		

	// construct the vector b
	Eigen::Matrix<fReal, nx*ny, 1> b;
	b.setZero();
	for(size_t i = 0; i < nx; ++i){
		for(size_t j = 0; j < ny; ++j){
			fReal uPlus, uMinus, vPlus, vMinus;
			uPlus = attributeTable["u"]->getValueAt(i + 1, j);
			uMinus = attributeTable["u"]->getValueAt(i, j);
			vPlus = attributeTable["v"]->getValueAt(i, j + 1);
			vMinus = attributeTable["v"]->getValueAt(i, j);
			b(j*nx + i) = -((uPlus - uMinus) * invGridLen + (vPlus - vMinus) * invGridLen);
		}
	}

	// density vector
	Eigen::Matrix<fReal, nx*ny, 1> rho;

	// solving Ax = b
	Eigen::MINRES<SparseMatrix<fReal>, Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner> minres;
    minres.compute(A);
    rho = minres.solve(b);

	// Populate updated quantities
    for(size_t i = 0; i < nx; ++i){
    	for(size_t j = 0; j < ny; ++j){
    		attributeTable["rho"]->setValueAt(i, j) = rho(j*nx + i);
    	}	
    }
    attributeTable["rho"]->swapBuffer();

    for(size_t i = 0; i < nx + 1; ++i){
    	for(size_t j = 0; j < ny; ++j){
    		if(i == 0){
    			attributeTable["u"]->setValueAt(i, j) = attributeTable["u"]->getValueAt(i, j) -
    			scale * invGridLen * (attributeTable["rho"]->getValueAt(i, j) - attributeTable["rho"]->getValueAt(nx - 1, j))
    		}
    		else if(i == nx){
    			attributeTable["u"]->setValueAt(i, j) = attributeTable["u"]->getValueAt(i, j) -
    			scale * invGridLen * (attributeTable["rho"]->getValueAt(0, j) - attributeTable["rho"]->getValueAt(i - 1, j))
    		}
    		else{
    			attributeTable["u"]->setValueAt(i, j) = attributeTable["u"]->getValueAt(i, j) -
    			scale * invGridLen * (attributeTable["rho"]->getValueAt(i, j) - attributeTable["rho"]->getValueAt(i - 1, j))
    		}
    	}
    }
    attributeTable["u"]->swapBuffer();

    for(size_t i = 0; i < nx; ++i){
    	for(size_t j = 0; j < ny + 1; ++j){
    		if(j == 0){
    			attributeTable["v"]->setValueAt(i, j) = attributeTable["v"]->getValueAt(i, j) -
    			scale * invGridLen * (attributeTable["rho"]->getValueAt(i, j) - attributeTable["v"]->getValueAt(i, ny - 1));
    		}
    		else if(j == ny){
    			attributeTable["v"]->setValueAt(i, j) = attributeTable["v"]->getValueAt(i, j) -
    			scale * invGridLen * (attributeTable["rho"]->getValueAt(i, j - 1) - attributeTable["v"]->getValueAt(i, 0));    			
    		}
    		else{
    			attributeTable["v"]->setValueAt(i, j) = attributeTable["v"]->getValueAt(i, j) -
    			scale * invGridLen * (attributeTable["rho"]->getValueAt(i, j) - attributeTable["v"]->getValueAt(i, j - 1));
    		}
    	}
    }
    attributeTable["v"]->swapBuffer();
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

fReal KaminoGrid::FBM(const fReal x, const fReal y){
    fReal total = 0.0f;
    fReal resolution = 10.f;
    fReal persistance = 0.5;
    int octaves = 4;

    for(int i = 0; i < octaves; i++){
        fReal freq = std::pow(2.0f, i);
        fReal amp = std::pow(persistance, i);
        total += amp * interpNoise2D(x * freq / resolution, y * freq / resolution);
    }
    fReal a = 1 - persistance;  // normalization

    return a * total / 2.0f;  // normalized, pseudorandom number between -1 and 1

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

void KaminoGrid::write_data_bgeo(const std::string& s, const int frame)
{
# ifndef _MSC_VER
    std::string file = s + std::to_string(frame) + ".bgeo";

    // TODO: interpolate velocities to grid centers and combine into vec2

    Partio::ParticlesDataMutable* parts = Partio::create();
    Partio::ParticleAttribute pH, vH;
    pH = parts->addAttribute("position", Partio::VECTOR, 3);
    vH = parts->addAttribute("v", Partio::VECTOR, 3);

    Eigen::Matrix<float, 3, 1> pos;
    Eigen::Matrix<float, 3, 1> vel;
    fReal velX, velY;

    for(size_t i = 0; i < nx; ++i){
        for(size_t j = 0; j < ny; ++j){
        	velX = (attributeTable["u"]->getValueAt(i, j) + attributeTable["u"]->getValueAt(i + 1, j)) / 2.0;
        	velY = (attributeTable["v"]->getValueAt(i, j) + attributeTable["v"]->getValueAt(i, j + 1)) / 2.0;
        	pos = Eigen::Matrix<float, 3, 1>(i * gridLen, j * gridLen, 0.0);
        	vel = Eigen::Matrix<float, 3, 1>(velX, velY, 0.0);
            int idx = parts->addParticle();
            float* p = parts->dataWrite<float>(pH, idx);
            float* v = parts->dataWrite<float>(vH, idx);
            for (int k = 0; k < 3; ++k){
                p[k] = pos(k, 0);
                v[k] = vel(k, 0);
            }
        }
    }
    Partio::write(file.c_str(), *parts);
    parts->release();
# endif
}

void KaminoGrid::initialize_velocity()
{
	// initialize u velocities
	fReal x = -gridLen / 2.0;
	fReal y = 0.0;
	fReal val = 0.0;

    for(size_t i = 0; i < nx; ++i){
        for(size_t j = 0; j < ny; ++j){
        	val = FBM(sin(2*M_PI*x / (nx*gridLen)), sin(2*M_PI*y / (ny*gridLen)));
            attributeTable["u"]->setValueAt(i, j, val);
            y += gridLen;
        }
        x += gridLen;
    }
    // initialize v velocities
    x = 0.0;
    y = -gridLen / 2.0;

    for(size_t i = 0; i < nx; ++i){
        for(size_t j = 0; j < ny + 1; ++j){
        	val = FBM(sin(2*M_PI*x / (nx*gridLen)), sin(2*M_PI*y / (ny*gridLen)));
            attributeTable["v"]->setValueAt(i, j, val);
            y += gridLen;
        }
        x += gridLen;
    }
}