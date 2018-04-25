# include "../include/KaminoQuantity.h"

KaminoParticles::KaminoParticles(int n, fReal particleDensity, fReal radius, KaminoSolver* solver) :
                            particleDensity(particleDensity), radius(radius), parentSolver(solver)
{
    fReal delta = 1 / particleDensity;
    fReal halfDelta = delta / 2.0;
    
    for(unsigned int i = 0; i < 2 * n; ++i){
        for(unsigned int j = 0; j < n; ++j){
            // distribute in phi and theta randomly
            // +/- is 50/50
            fReal signPhi = static_cast <fReal> (rand()) / static_cast <fReal> (RAND_MAX);
            signPhi = signPhi >= 0.5 ? 1.0 : -1.0;
            fReal signTheta = static_cast <fReal> (rand()) / static_cast <fReal> (RAND_MAX);
            signTheta = signTheta >= 0.5 ? 1.0 : -1.0;

            // get random value between 0 and halfDelta in +/- direction
            fReal randPhi = signPhi * halfDelta * static_cast <fReal> (rand()) / static_cast <fReal> (RAND_MAX);
            fReal randTheta = signTheta * halfDelta * static_cast <fReal> (rand()) / static_cast <fReal> (RAND_MAX);

            // assign positions (phi, theta)
            fReal phi = i * delta + randPhi;
            fReal theta = j * delta + randTheta;
            Eigen::Matrix<fReal, 2, 1> pos(phi, theta);
            positions.push_back(pos);
        }
    }
}

KaminoParticles::~KaminoParticles()
{
}

void KaminoParticles::updatePositions(KaminoQuantity* u, KaminoQuantity* v, fReal deltaT)
{
# ifdef OMParallelize
# pragma omp parallel for
# endif
    for(unsigned int i = 0; i < positions.size(); ++i){
        fReal uPhi = u->sampleAt(positions[i][0], positions[i][1], parentSolver->uNorthP, parentSolver->uSouthP);
        fReal uTheta = v->sampleAt(positions[i][0], positions[i][1], parentSolver->uNorthP, parentSolver->uSouthP);
        fReal nextPhi;
        fReal nextTheta;
        // problem at the south pole
        if(positions[i][1] > M_PI - 1E-7 && positions[i][1] < M_PI + 1E-7){
            nextPhi = positions[i][0] + uPhi * deltaT / (radius * sin(1E-10));
        }
        // problem at the north pole
        else if(positions[i][1] < 1E-7 && positions[i][1] > -1E-7){
            nextPhi = positions[i][0] + uPhi * deltaT / (radius * sin(1E-10));
        }
        else{
            nextPhi = positions[i][0] + uPhi * deltaT / (radius * sin(positions[i][1]));
        } 
        nextTheta = positions[i][1] + uTheta * deltaT / radius;
        // wrap particles
        bool check = validatePhiTheta(nextPhi, nextTheta);
        positions[i][0] = nextPhi;
        positions[i][1] = nextTheta;
    }
}

void KaminoParticles::write_data_bgeo(const std::string& s, const int frame)
{
//# ifndef _MSC_VER
    std::string file = s + std::to_string(frame) + ".bgeo";
    //std::cout << "Writing to: " << file << std::endl;
    Partio::ParticlesDataMutable* parts = Partio::create();
    Partio::ParticleAttribute posH, vH;
    posH = parts->addAttribute("position", Partio::VECTOR, 3);
    vH = parts->addAttribute("v", Partio::VECTOR, 3);

    for(unsigned int i = 0; i < positions.size(); ++i){
        int idx = parts->addParticle();
        float* p = parts->dataWrite<float>(posH, idx);
        float* v = parts->dataWrite<float>(vH, idx);
        Eigen::Matrix<float, 3, 1> pos(positions[i][0], positions[i][1], 0.0);
        mapPToSphere(pos);
        for (int k = 0; k < 3; ++k){
            p[k] = pos(k, 0);
            v[k] = 0.0;
        }
    }
    Partio::write(file.c_str(), *parts);
    parts->release();
//# endif
}

void KaminoParticles::mapPToSphere(Eigen::Matrix<float, 3, 1>& pos) const
{
    float theta = pos[1];
    float phi = pos[0];
    pos[0] = radius * sin(theta) * cos(phi);
    pos[2] = radius * sin(theta) * sin(phi);
    pos[1] = radius * cos(theta);
}