# include "../include/KaminoQuantity.h"

KaminoSolver::KaminoSolver(size_t nx, size_t ny, fReal gridLength, fReal frameDuration) :
	nx(nx), ny(ny), gridLen(gridLength), frameDuration(frameDuration),
	timeStep(0.0), timeElapsed(0.0)
{
	addAttr("u", 0.5, 0.0);		// u velocity
	addAttr("v", 0.0, 0.5);		// v velocity
	addAttr("p");				// p density

	initialize_velocity();
}

KaminoSolver::~KaminoSolver()
{
	for (auto& attr : this->attributeTable)
	{
		delete attr.second;
	}
}

void KaminoSolver::stepForward(fReal timeStep)
{
	this->timeStep = timeStep;
	advection();
	// bodyForce();
	projection();
# ifdef DEBUGBUILD
	/*for (unsigned gridX = 0; gridX != nx; ++gridX)
	{
		for (unsigned gridY = 0; gridY != ny; ++gridY)
		{
			std::cout << attributeTable["u"]->getValueAt(gridX, gridY) << '\t';
		}
		std::cout << '\n';
	}*/
# endif
}

void KaminoSolver::addAttr(std::string name, fReal xOffset, fReal yOffset)
{
	KaminoQuantity* ptr = new KaminoQuantity(name, this->nx, this->ny, this->gridLen, xOffset, yOffset);
	this->attributeTable.emplace(std::pair<std::string, KaminoQuantity*>(name, ptr));
}

KaminoQuantity* KaminoSolver::getAttributeNamed(std::string name)
{
	return (*this)[name];
}

KaminoQuantity* KaminoSolver::operator[](std::string name)
{
	return this->attributeTable.at(name);
}

fReal KaminoSolver::FBM(const fReal x, const fReal y) {
	fReal total = 0.0f;
	fReal resolution = 1.0;
	fReal persistance = 0.5;
	int octaves = 4;

	for (int i = 0; i < octaves; i++) {
		fReal freq = std::pow(2.0f, i);
		fReal amp = std::pow(persistance, i);
		total += amp * interpNoise2D(x * freq / resolution, y * freq / resolution);
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

void KaminoSolver::write_data_bgeo(const std::string& s, const int frame)
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

	for (size_t i = 0; i < nx; ++i) {
		for (size_t j = 0; j < ny; ++j) {
			velX = (attributeTable["u"]->getValueAt(i, j) + attributeTable["u"]->getValueAt(i + 1, j)) / 2.0;
			velY = (attributeTable["v"]->getValueAt(i, j) + attributeTable["v"]->getValueAt(i, j + 1)) / 2.0;
			pos = Eigen::Matrix<float, 3, 1>(i * gridLen, j * gridLen, 0.0);
			vel = Eigen::Matrix<float, 3, 1>(velX, velY, 0.0);
			int idx = parts->addParticle();
			float* p = parts->dataWrite<float>(pH, idx);
			float* v = parts->dataWrite<float>(vH, idx);
			for (int k = 0; k < 3; ++k) {
				p[k] = pos(k, 0);
				v[k] = vel(k, 0);
			}
		}
	}
	Partio::write(file.c_str(), *parts);
	parts->release();
# endif
}

void KaminoSolver::initialize_velocity()
{
	// initialize u velocities
	fReal x = -gridLen / 2.0;
	fReal y = 0.0;
	fReal val = 0.0;

	for (size_t i = 0; i < nx; ++i) {
		for (size_t j = 0; j < ny; ++j) {
			val = FBM(sin(2 * M_PI*x / (nx*gridLen)), sin(2 * M_PI*y / (ny*gridLen)));
			attributeTable["u"]->setValueAt(i, j, val);
			y += gridLen;
		}
		x += gridLen;
	}
	// initialize v velocities
	x = 0.0;
	y = -gridLen / 2.0;

	for (size_t i = 0; i < nx; ++i) {
		for (size_t j = 0; j < ny; ++j) {
			val = FBM(std::sin(2.0 * M_PI *x / (nx * gridLen)), sin(2.0 * M_PI * y / (ny * gridLen)));
			attributeTable["v"]->setValueAt(i, j, val);
			y += gridLen;
		}
		x += gridLen;
	}
}

void KaminoSolver::advection()
{
	for (auto quantity : this->attributeTable)
	{
		KaminoQuantity* attr = quantity.second;
		for (size_t gridX = 0; gridX < this->nx; ++gridX)
		{
			for (size_t gridY = 0; gridY < this->ny; ++gridY)
			{
				fReal gX = attr->getXCoordinateAt(gridX);
				fReal gY = attr->getYCoordinateAt(gridY);

				fReal uG = (*this)["u"]->sampleAt(gX, gY);
				fReal vG = (*this)["v"]->sampleAt(gX, gY);

				fReal midX = gX - 0.5 * timeStep * uG;
				fReal midY = gY - 0.5 * timeStep * vG;

				fReal uMid = (*this)["u"]->sampleAt(midX, midY);
				fReal vMid = (*this)["v"]->sampleAt(midX, midY);

				fReal pX = gX - timeStep * uMid;
				fReal pY = gY - timeStep * vMid;

				fReal advectedVal = attr->sampleAt(pX, pY);
				attr->writeValueTo(gridX, gridY, advectedVal);
			}
		}
	}
	this->swapAttrBuffers();
}

void KaminoSolver::projection()
{
	fReal density = 1.0;	// rest fluid density
	fReal scale = timeStep / density;
	fReal invGridLen = 1 / gridLen;

	// construct the matrix A
	// present construction assumes 2D fluid in every cell and toroidal BCs
	//Eigen::MatrixXf A(nx*ny, nx*ny);
	Eigen::SparseMatrix<fReal> A(nx*ny, nx*ny);
	A.setZero();

	// construct A row-by-row
	size_t k = 0;
	Eigen::VectorXd ARow(nx * ny);
	//Eigen::VectorXf ARow(nx * ny);
	for(size_t i = 0; i < nx; ++i){
		for(size_t j = 0; j < ny; ++j){
			ARow.setZero();
			ARow(j*nx + i) = 4;
			i > (nx - 2) ? (ARow(j*nx) = -1) : (ARow(j*nx + i + 1) = -1);
			i < 1 ? (ARow(j*nx + nx - 1) = -1) : (ARow(j*nx + i - 1) = -1);
			j > (ny - 2) ? (ARow(i) = -1) : (ARow((j + 1)*nx + i) = -1);
			j < 1 ? (ARow((ny - 1)*nx + i) = -1) : (ARow((j - 1)*nx + i) = -1);
			for(int l = 0; l < nx * ny; ++l){
				A.coeffRef(k, l) = ARow(l);
			}
			//A.row(k) = ARow;
			k++;
		}
	}
	A *= scale;		

	// construct the vector b
	Eigen::VectorXd b(nx * ny);
	//Eigen::VectorXf b(nx * ny);
	b.setZero();
	for(size_t i = 0; i < nx; ++i){
		for(size_t j = 0; j < ny; ++j){
			fReal uPlus, uMinus, vPlus, vMinus;
			i > (nx - 2) ? (uPlus = attributeTable["u"]->getValueAt(0, j)) : (uPlus = attributeTable["u"]->getValueAt(i + 1, j));
			uMinus = attributeTable["u"]->getValueAt(i, j); 
			j > (ny - 2) ? (vPlus = attributeTable["v"]->getValueAt(i, 0)) : (vPlus = attributeTable["v"]->getValueAt(i, j + 1));
			vMinus = attributeTable["v"]->getValueAt(i, j);
			b(j*nx + i) = -((uPlus - uMinus) * invGridLen + (vPlus - vMinus) * invGridLen);
		}
	}

	// pressure vector
	//Eigen::VectorXd p(nx * ny);
	//Eigen::VectorXf p(nx * ny);
	//p.setZero();

	// solving Ax = b
	// Eigen::MINRES<Eigen::MatrixXf, Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner> minres;
 	// minres.compute(A);
 	// p = minres.solve(b);

	Eigen::ConjugateGradient<Eigen::SparseMatrix<fReal>, Eigen::Lower|Eigen::Upper> cg(A);
	//cg.compute(A);
	Eigen::VectorXd p = cg.solve(b);

	// Populate updated pressure values
    for(size_t i = 0; i < nx; ++i){
    	for(size_t j = 0; j < ny; ++j){
    		attributeTable["p"]->writeValueTo(i, j, p(j*nx + i));
    	}	
    }

    // Populate updated u values
    for(size_t i = 0; i < nx; ++i){
    	for(size_t j = 0; j < ny; ++j){
    		if(i == 0){
    			attributeTable["u"]->writeValueTo(i, j, attributeTable["u"]->getValueAt(i, j) -
    			scale * invGridLen * (attributeTable["p"]->getNextValueAt(i, j) - attributeTable["p"]->getNextValueAt(nx - 1, j)));
    		}
    		else{
    			attributeTable["u"]->writeValueTo(i, j, attributeTable["u"]->getValueAt(i, j) -
    			scale * invGridLen * (attributeTable["p"]->getNextValueAt(i, j) - attributeTable["p"]->getNextValueAt(i - 1, j)));
    		}
    	}
    }

    // Populate updated v values
    for(size_t i = 0; i < nx; ++i){
    	for(size_t j = 0; j < ny; ++j){
    		if(j == 0){
    			attributeTable["v"]->writeValueTo(i, j, attributeTable["v"]->getValueAt(i, j) -
    			scale * invGridLen * (attributeTable["p"]->getNextValueAt(i, j) - attributeTable["v"]->getNextValueAt(i, ny - 1)));
    		}
    		else{
    			attributeTable["v"]->writeValueTo(i, j, attributeTable["v"]->getValueAt(i, j) -
    			scale * invGridLen * (attributeTable["p"]->getNextValueAt(i, j) - attributeTable["v"]->getNextValueAt(i, j - 1)));
    		}
    	}
    }
    this->swapAttrBuffers();
}

void KaminoSolver::swapAttrBuffers()
{
	for (auto quantity : this->attributeTable)
	{
		quantity.second->swapBuffer();
	}
}