# include "../include/KaminoQuantity.h"

int main(int argc, char** argv)
{
	KaminoSolver solver = KaminoSolver(2, 2, 0.1);
	for (unsigned i = 0; i != 100; ++i)
	{
		solver.stepForward(0.001);
		
	}
}