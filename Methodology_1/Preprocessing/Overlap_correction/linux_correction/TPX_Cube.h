#pragma once

#include "FitsFile.h"

#define MAX_SHUTTERS_AVAILABLE 256

class TPX_Cube
{
public:
	TPX_Cube(char *Filename_par);
	TPX_Cube(void);
	~TPX_Cube(void);

	char			*FileNameBase;	//file name before the indexes
	char			*Dir_Name;	//file name before the indexes
	int				nSlices;		//number of spectral elements
	unsigned short	nShtrs;			//number of shutters used for acquisition
	unsigned long	ShtrCounts[MAX_SHUTTERS_AVAILABLE];//array of count values in each of 255 chutters available
	double			ShtrStart[MAX_SHUTTERS_AVAILABLE]; //array of start times for shutters
	double			ShtrEnd[MAX_SHUTTERS_AVAILABLE];   //array of end times for shutters

	double			*TOF_Arr;		//array of TOF values
	unsigned short	*Shtr_Arr;		//array of shutter number for each slice
	double			*TimeBin_Arr;	//array of time bin values for each slice
	double			*Wvl_Arr;		//array of wavelength values, angstroms
	double			*E_Arr;			//array of neutron energies, eV
	double			L;				//flight path, initially is set to 0, that means no Wvl_Arr is calculated yet
	double			Tshift;			//s, shift of trigger used in the experiment

	double			*Spectra_Arr;	//array of spectra values calculated last time, at Constructor - total N of counts in image
									//the array is created at the same time as TOF, same dimension
	unsigned short	SpectrX0;		//current postion of the spectra area calculated
	unsigned short	SpectrX1;		//current postion of the spectra area calculated
	unsigned short	SpectrY0;		//current postion of the spectra area calculated
	unsigned short	SpectrY1;		//current postion of the spectra area calculated
	int				SpectrIndLeft;	//left index of spectra calculated 
	int				SpectrIndRight;	//Right index of spectra calculated 

	FitsFile		*Slices;		//images in spectral cube

	int		LoadFromFiles(char *FileName_par);
	int		SaveCubeToFitsFiles(char *FileName_par);
	int		CreateCopy(int ImgTypeAST, TPX_Cube *SrcCube); //create copy of existing Cube, but the image file type can be changed - need that for the floating point of normalized of overlap corrected cube
	int		SaveOverlapCorrectedCube(char *FileNameBase_par, bool ROTAX);
	double	*GetSpectra(unsigned short *iDim, unsigned short X0, unsigned short X1, unsigned short Y0, unsigned short Y1, double TOF0, double TOF1);
	double 	TOF_To_Wvlngth(double TOF_par, double FlightPath, double Tshift_par); //converts TOF to Avlngth
	double 	Wvlngth_to_TOF(double Wvlgth, double FlightPath, double Tshift_par);  //converts wavelength to TOF measured (including the Tshift value)
	double	TOF_To_E(double TOF_par, double FlightPath_par, double Tshift_par);
	double	E_to_TOF(double E, double FlightPath_par, double Tshift_par);
	int		SetExperimVariables(double FlightPath_par, double Tshift_par);
	int		SaveLastSpectra(char *FileName);
	double	GetNearestArrValue(double Val, double *Arr, int iDim);


//	inline int round(double X);

//protected:
//    void initHwInfoItems();

//private:
//	ofstream	_fLog;						// Output stream for log messages AST};
};