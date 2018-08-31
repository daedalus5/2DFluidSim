# include "../include/HH16Quantity.h"
# include "../include/CubicSolver.h"
# include "../include/KaminoTimer.h"

// CONSTRUCTOR / DESTRUCTOR >>>>>>>>>>

HH16Solver::HH16Solver(size_t nPhi, size_t nTheta, fReal radius, fReal gridLength, fReal frameDuration) :
    nPhi(nPhi), nTheta(nTheta), radius(radius), gridLen(gridLength), invGridLen(1.0 / gridLength), frameDuration(frameDuration),
    timeStep(0.0), timeElapsed(0.0), advectionTime(0.0f), geometricTime(0.0f), projectionTime(0.0f)
{
    addAttr("u");         // u velocity
    addAttr("v");         // v velocity
    addAttr("p");         // p pressure
    addAttr("density");   // density
    
    initialize_velocity();
    initialize_pressure();
    initialize_density();
}

HH16Solver::~HH16Solver()
{
    for (auto& attr : this->attr)
    {
        delete attr.second;
    }

    float totalTimeUsed = this->advectionTime + this->geometricTime + this->projectionTime;
    std::cout << "Total time used for advection : " << this->advectionTime << std::endl;
    std::cout << "Total time used for geometric : " << this->geometricTime << std::endl;
    std::cout << "Total time used for projection : " << this->projectionTime << std::endl;
    std::cout << "Percentage of advection : " << advectionTime / totalTimeUsed * 100.0f << "%" << std::endl;
    std::cout << "Percentage of geometric : " << geometricTime / totalTimeUsed * 100.0f << "%" << std::endl;
    std::cout << "Percentage of projection : " << projectionTime / totalTimeUsed * 100.0f << "%" << std::endl;
}

void HH16Solver::stepForward(fReal timeStep)
{
    this->timeStep = timeStep;
    //advectionScalar();
    KaminoTimer timer;
    timer.startTimer();
    advectionSpeed();
    this->advectionTime += timer.stopTimer();
    
    this->swapAttrBuffers();

    timer.startTimer();
    geometric(); // Buffer is swapped here
    this->geometricTime += timer.stopTimer();
    //bodyForce();
    timer.startTimer();
    projection();
    this->projectionTime += timer.stopTimer();
    this->timeElapsed += timeStep;
}

// In KaminoSolver.cpp
/*
// Phi: 0 - 2pi  Theta: 0 - pi
bool validatePhiTheta(fReal & phi, fReal & theta)
{
    int loops = static_cast<int>(std::floor(theta / M_2PI));
    theta = theta - loops * M_2PI;
    // Now theta is in 0-2pi range

    bool isFlipped = false;
    
    if (theta > M_PI)
    {
        theta = M_2PI - theta;
        phi += M_PI;
        isFlipped = true;
    }

    loops = static_cast<int>(std::floor(phi / M_2PI));
    phi = phi - loops * M_2PI;
    // Now phi is in 0-2pi range

    return isFlipped;
}
*/

void HH16Solver::bodyForce()
{
    fReal gravity = 9.8;
    HH16Quantity* v = attr["v"];

    for(size_t j = 0; j < nTheta + 1; ++j){
        for(size_t i = 0; i < nPhi; ++i){
            fReal vBeforeUpdate = v->getValueAt(i, j);
            fReal theta = j*gridLen;
            v->writeValueTo(i, j, vBeforeUpdate + gravity * sin(theta) * timeStep);
        }
    }

    v->swapBuffer();
}

/* Duplicate of getIndex() in HH16Quantity */
size_t HH16Solver::getIndex(size_t x, size_t y)
{
    return y * nPhi + x;
}

void HH16Solver::addAttr(std::string name)
{
    size_t attrnPhi = this->nPhi;
    size_t attrnTheta = this->nTheta;

    HH16Quantity* ptr = new HH16Quantity(name, attrnPhi, attrnTheta, this->gridLen);
    this->attr.emplace(std::pair<std::string, HH16Quantity*>(name, ptr));
}

HH16Quantity* HH16Solver::getAttributeNamed(std::string name)
{
    return (*this)[name];
}

HH16Quantity* HH16Solver::operator[](std::string name)
{
    if (this->attr.find(name) == this->attr.end())
    {
        std::cerr << "attribute not found in map" << std::endl;
    }
    return this->attr.at(name);
}

void HH16Solver::swapAttrBuffers()
{
    for (auto quantity : this->attr)
    {
        quantity.second->swapBuffer();
    }
}


// <<<<<<<<<<
// OUTPUT >>>>>>>>>>


void HH16Solver::write_data_bgeo(const std::string& s, const int frame)
{
    std::string file = s + std::to_string(frame) + ".bgeo";
    std::cout << "Writing to: " << file << std::endl;

    Partio::ParticlesDataMutable* parts = Partio::create();
    Partio::ParticleAttribute pH, vH, psH, dens;
    pH = parts->addAttribute("position", Partio::VECTOR, 3);
    vH = parts->addAttribute("v", Partio::VECTOR, 3);
    psH = parts->addAttribute("pressure", Partio::VECTOR, 1);
    dens = parts->addAttribute("density", Partio::VECTOR, 1);

    Eigen::Matrix<float, 3, 1> pos;
    Eigen::Matrix<float, 3, 1> vel;
    fReal pressure, densityValue;
    fReal velX, velY;

    HH16Quantity* u = attr["u"];
    HH16Quantity* v = attr["v"];
    fReal uRight, uLeft, vUp, vDown;

    size_t upi, vpi;

    for (size_t j = 0; j < nTheta; ++j) {
        for (size_t i = 0; i < nPhi; ++i) {
            uLeft = u->getValueAt(i, j);
            i == (nPhi - 1) ? upi = 0 : upi = i + 1;
            vDown = v->getValueAt(i, j);
            j == (nTheta - 1) ? vpi = 0 : vpi = j + 1;
            uRight = u->getValueAt(upi, j);
            vUp = u->getValueAt(i, vpi);

            velX = (uLeft + uRight) / 2.0;
            velY = (vUp + vDown) / 2.0;

            pos = Eigen::Matrix<float, 3, 1>(i * gridLen, j * gridLen, 0.0);
            vel = Eigen::Matrix<float, 3, 1>(0.0, velY, velX);
            mapVToSphere(pos, vel);
            mapPToSphere(pos);

            pressure = attr["p"]->getValueAt(i, j);
            densityValue = attr["density"]->getValueAt(i, j);
            
            int idx = parts->addParticle();
            float* p = parts->dataWrite<float>(pH, idx);
            float* v = parts->dataWrite<float>(vH, idx);
            float* ps = parts->dataWrite<float>(psH, idx);
            float* de = parts->dataWrite<float>(dens, idx);

            ps[0] = 1000.0 * radius * pressure / timeStep;
            de[0] = densityValue;

            for (int k = 0; k < 3; ++k) {
                p[k] = pos(k, 0);
                v[k] = vel(k, 0);
            }
        }
    }

    Partio::write(file.c_str(), *parts);
    parts->release();
}

void HH16Solver::mapPToSphere(Eigen::Matrix<float, 3, 1>& pos) const
{
    float theta = pos[1];
    float phi = pos[0];
    pos[0] = radius * sin(theta) * cos(phi);
    pos[2] = radius * sin(theta) * sin(phi);
    pos[1] = radius * cos(theta);
}

void HH16Solver::mapVToSphere(Eigen::Matrix<float, 3, 1>& pos, Eigen::Matrix<float, 3, 1>& vel) const
{
    float theta = pos[1];
    float phi = pos[0];

    float u_theta = vel[1];
    float u_phi = vel[2];

    vel[0] = cos(theta) * cos(phi) * u_theta - sin(phi) * u_phi;
    vel[2] = cos(theta) * sin(phi) * u_theta + cos(phi) * u_phi;
    vel[1] = -sin(theta) * u_theta;
}

void HH16Solver::mapToCylinder(Eigen::Matrix<float, 3, 1>& pos) const
{
    //float radius = 5.0;
    float phi = 2*M_PI*pos[0] / (nPhi * gridLen);
    float z = pos[1];
    pos[0] = radius * cos(phi);
    pos[1] = radius * sin(phi);
    pos[2] = z;
}