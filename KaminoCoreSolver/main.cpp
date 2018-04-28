# include "include/Kamino.h"
# include <fstream>

int main(int argc, char** argv)
{
	std::string configFile = argv[1];
	std::fstream fin;
	fin.open(configFile, std::ios::in);
	fReal radius; size_t nTheta; fReal particleDensity;
	float dt; float DT; int frames;
	std::string gridPath; std::string particlePath;
	std::string densityImage; std::string solidImage; std::string colorImage;

	fin >> radius;
	fin >> nTheta;
	fin >> particleDensity;
	fin >> dt;
	fin >> DT;
	fin >> frames;
	fin >> gridPath;
	fin >> particlePath;

	fin >> densityImage;
	if (densityImage == "null")
	{
		densityImage = "";
	}
	fin >> solidImage;
	if (solidImage == "null")
	{
		solidImage == "";
	}
	fin >> colorImage;
	if (colorImage == "null")
	{
		colorImage = "";
	}

    Kamino KaminoInstance(radius, nTheta, particleDensity, dt, DT, frames,
		gridPath, particlePath, densityImage, solidImage, colorImage);
    KaminoInstance.run();
    return 0;
}