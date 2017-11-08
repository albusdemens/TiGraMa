#include <string.h>
#include "TPX_Cube.h"
#include "../common/utils.h"

using namespace std;


int main(int argc, char *argv[])
{
	printf("Overlap correction code. \nProgram creates a new directory with the corrected data\n");
	if (argc < 2)
	{
		printf("usage: program [fits file] (/ROTAX)\n");
		return 0;
	}

	bool ROTAX = false;
	if (argc > 2 && strncmp(argv[2], "/ROTAX", MAX(strlen(argv[2]), strlen("/ROTAX"))) == 0)
	{
		ROTAX = true;
		printf("ROTAX-specific correction will be performed!!\n");
		printf("Will correct the shutter values by 4/5 as ROTAX trigger comes even if the pulse is sent to TS2\n");
	}

	TPX_Cube *Cube = new TPX_Cube(argv[1]);
	Cube->SaveOverlapCorrectedCube(Cube->FileNameBase, ROTAX);
	unsigned short iTmp;
	double* dTmp = Cube->GetSpectra(&iTmp, 10, 50, 10, 50, 5e-3, 10e-3);
	return 0;
}

