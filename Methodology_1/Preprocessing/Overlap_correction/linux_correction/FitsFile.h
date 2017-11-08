#pragma once

#include <fstream>    
#include <sstream>	  
#include <iostream>	  

#include "/usr/local/Cellar/cfitsio/3.370/include/fitsio.h"  /* required by every program that uses CFITSIO  */


class FitsFile
{
public:
	FitsFile(void);
    ~FitsFile();

	char			*FileName;		//file name of the current file in memory
	unsigned short	FileIndex;		//Index within the filename
	int				iX;				//dimension of X axis
	int				iY;				//dimension of Y axis
	char			*Img;			//fits image, need to convert it to type once it is known
	int				ImgTypeAST;		//defines the type of image: 0 - byte, 1  short, 2 - ushort, 3 - int, 4 - long, 5 - ulong, 6 - float, 7- double

	float			TOF;			//time of flight, in seconds
	float			TimeBin;		//width of the time bin, in seconds
	unsigned long	Counts;			//counts in the image
	unsigned int	Triggers;		//number of triggers received by that image
	unsigned short	ShtrN;			//shutter number at which the image was acquired, i.e. TOF belongs to shtrN


	int		CreateEmptyImage(int ImgType_par, int iX_par, int iY_par);
	int		CopyImg(FitsFile *ImgToCopy);
	int		AddImg(FitsFile *ImgToAdd);
	int		DivideByImage(FitsFile *ImgToAdd);
	int		ZeroImage(void);
	int		DivideByConst(double Divider);
	double 	getCountsInArea(unsigned short X0, unsigned short X1, unsigned short Y0, unsigned short Y1);

	int		ReadFitsFile(char *FileName_par);
	int		SaveImage_to_FITS(char *FileName_par);
	int		getImgBytesPerPixel(int ImgType_par);

//	inline int round(double X);

//protected:
//    void initHwInfoItems();

//private:
//	ofstream	_fLog;						// Output stream for log messages AST
};

enum FitsErrorCodes
{
	FITS_WRITE_ERROR		= -0x10000001,	//error wrting FITS file
	FITS_READ_ERROR			= -0x10000002,	//error reading FITS file
	FITS_COPYIMG_ERROR		= -0x10000003,	//error copying image
	FITS_ADDIMG_ERROR		= -0x10000004,	//error adding image
	FITS_DIVIMG_ERROR		= -0x10000005,	//error dividing image
	CUBE_CONSTRUCT_ERR		= -0x10000006,	//error constructing TPX cube
	CUBE_N_SLICES_ERR		= -0x10000007,	//number of slices in the cube cannot be <=0
	CUBE_OVERLAP_CORR_ERR	= -0x10000008,	//error in correcting overlaps
	CUBE_FLIGHT_PATH_ERR	= -0x10000009,	//flight path cannot be 0 or negative
	CUBE_BUILD_SPECTRA_ERR	= -0x10000010,	//building spectra called with wrong parameters
	FITS_BUILD_SPECTRA_ERR	= -0x10000011,	//building spectra called with wrong parameters
};
