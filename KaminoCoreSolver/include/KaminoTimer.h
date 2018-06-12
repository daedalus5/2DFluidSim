# pragma once

# include <chrono>

class KaminoTimer
{
private:
	std::chrono::high_resolution_clock::time_point tStart;
	std::chrono::high_resolution_clock::time_point tStop;
public:
	void startTimer();
	float stopTimer();
};