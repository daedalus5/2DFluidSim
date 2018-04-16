# include "../include/KaminoQuantity.h"

KaminoParticles::KaminoParticles(int n, fReal radius) : radius(radius)
{
    fReal delta = M_PI / n;
    Eigen::Matrix<fReal, 2, 1> pos(0.0, 0.0);
    for(unsigned int i = 0; i < 2 * n; ++i){
        for(unsigned int j = 0; j < n; ++j){
            pos(0, 0) = i * delta;
            pos(1, 0) = j * delta;
            positions.push_back(pos);
        }
    }
}

KaminoParticles::~KaminoParticles()
{
}

void KaminoParticles::updatePositions(KaminoQuantity* u, KaminoQuantity* v, fReal deltaT)
{
    for(unsigned int i = 0; i < positions.size(); ++i){
        fReal uPhi = u->sampleAt(positions[i][0], positions[i][1]);
        fReal uTheta = v->sampleAt(positions[i][0], positions[i][1]);
        positions[i][0] += uPhi * deltaT / (radius * sin(positions[i][1]));
        positions[i][1] += uTheta * deltaT / radius;
    }
}

void KaminoParticles::write_data_bgeo(const std::string& s, const int frame)
{
    std::string file = s + std::to_string(frame) + ".bgeo";
    Partio::ParticlesDataMutable* parts = Partio::create();
    Partio::ParticleAttribute posH, vH;
    posH = parts->addAttribute("position", Partio::VECTOR, 3);
    vH = parts->addAttribute("v", Partio::VECTOR, 3);

    for(unsigned int i = 0; i < positions.size(); ++i){
        int idx = parts->addParticle();
        float* p = parts->dataWrite<float>(posH, idx);
        float* v = parts->dataWrite<float>(vH, idx);
        for (int k = 0; k < 3; ++k){
            p[k] = positions[i][k];
            v[k] = 0.0;
        }
    }
    Partio::write(file.c_str(), *parts);
    parts->release();
}