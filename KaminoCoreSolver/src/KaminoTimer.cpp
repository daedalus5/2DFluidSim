# include "../include/KaminoTimer.h"

void KaminoTimer::startTimer()
{
	this->tStart = std::chrono::high_resolution_clock::now();
}

float KaminoTimer::stopTimer()
{
	this->tStop = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> timeElapsed = 
		std::chrono::duration_cast<std::chrono::duration<float>>(tStop - tStart);
	return timeElapsed.count();
}