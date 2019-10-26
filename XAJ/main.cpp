#include "XAJModel.h"

void main() {
	XAJ xaj;
	xaj.Initialize(24, 537, 2, 1, "input.txt", 0, 70, 80, 0.1, 20, 40, 20);
	float Parameters[14] = { 0.65, 0.11, 0, 20, 75, 80, 3, 20, 1, 0.3, 0.41, 0.99, 0.6, 3 };
	float UH[3] = {0.3, 0.6, 0.1};
	xaj.SetParameters(Parameters, UH);
	xaj.RunModel();
	xaj.Routing();
	xaj.SaveOutput("output.csv");
}