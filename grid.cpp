#include "grid.h"

Grid::Grid(int nx, int ny, float dx, float dy) : nx(nx), ny(ny), dx(dx), dy(dy)
{
    u = new float[nx][ny];
    v = new float[nx][ny];
    positions = new Eigen::Matrix<float, 2, 1>[nx][ny];
    velocities = new Eigen::Matrix<float, 2, 1>[nx][ny];
    build_particle_grid();
    distribute_velocity();
}

Grid::~Grid(){}

void Grid::build_particle_grid()
{
    float x = 0.f;
    float y = 0.f;

    for(int i = 0; i < ny; ++i){
        for(int j = 0; j < nx; ++j){
            positions[i][j] = Eigen::Matrix<float, 2, 1>(x, y);
            x += dx;
        }
        y += dy;
    }
}

void Grid::distribute_velocity()
{
    for(int i = 0; i < nx; ++i){
        for(int j = 0; j < ny; ++j){
            velocities[i][j] = Eigen::Matrix<float, 2, 1>(0.0, 0.0);
        }
    }
}

void Grid::write_data_bgeo(const std::string& s, const int frame)
{
    std::string file = s + std::to_string(frame) + ".bgeo";
    Partio::ParticlesDataMutable* parts = Partio::create();
    Partio::ParticleAttribute pH, vH;
    pH = parts->addAttribute("p", Partio::VECTOR, 2);
    vH = parts->addAttribute("v", Partio::VECTOR, 2);

    for(int i = 0; i < ny; ++i){
        for(int j = 0; j < nx; ++j){
            int idx = parts->addParticle();
            float* p = parts->dataWrite<float>(pH, idx);
            float* v = parts->dataWrite<float>(vH, idx);
            for (int k = 0; k < 2; ++k){
                p[k] = positions[i][j][k];
                v[k] = velocities[i][j][k];
            }
        }
    }
    Partio::write(file.c_str(), *parts);
    parts->release();
}

Eigen::Matrix<float, 2, 1> Grid::FBM() const
{
    return Eigen::Matrix<float, 2, 1>(0.0,0.0);
}