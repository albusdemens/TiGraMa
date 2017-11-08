//#include "StdAfx.h"
#include <string.h>
#include <memory.h>
#include "FitsFile.h"

using namespace std;  

#define sprintf_s(a, b, c, d) snprintf(a, b, c, d)
#define memcpy_s(a, b, c, d) memcpy(a, c, d)
#define sscanf_s sscanf


#define DEBUG_LOG_CONSOLE
//#define DEBUG_LOG_FILE

#ifdef DEBUG_LOG_FILE
	#define LOGCALL(a)			_fLog<<a<<endl
	#define LOGCALL1(a,b)		_fLog<<a <<", "<<b<<endl
	#define LOGCALL2(a,b,c)		_fLog<<a <<", "<<b<<", "<<c<<endl
	#define LOGCALL3(a,b,c,d)	_fLog<<a <<", "<<b<<", "<<c<<", "<<d<<endl
	#define LOGCALL4(a,b,c,d,e)	_fLog<<a <<", "<<b<<", "<<c<<", "<<d<<","<<e<<endl
#else 
	#ifdef DEBUG_LOG_CONSOLE
		#define LOGCALL(a)			cout<<a<<endl
		#define LOGCALL1(a,b)		cout<<a <<", "<<b<<endl
		#define LOGCALL2(a,b,c)		cout<<a <<", "<<b<<", "<<c<<endl
		#define LOGCALL3(a,b,c,d)	cout<<a <<", "<<b<<", "<<c<<", "<<d<<endl
		#define LOGCALL4(a,b,c,d,e) cout<<a <<", "<<b<<", "<<c<<", "<<d<<","<<e<<endl
	#else
		#define LOGCALL(a)
		#define LOGCALL1(a,b)
		#define LOGCALL2(a,b,c)
		#define LOGCALL3(a,b,c,d)
	#endif
#endif

/*--------------------------------------------------------------------*/
FitsFile::FitsFile(void)
{
	#ifdef DEBUG_LOG_FILE
	_fLog.open("FitsFile.txt");
	#endif
	//LOGCALL("Initialize FitsFile");

	FileName=NULL;	//file name of the current file in memory
	FileIndex=0;	//Index within the filename
	iX=0;			//dimension of X axis
	iY=0;			//dimension of Y axis
	Img=NULL;		//fits image, need to convert it to type once it is known
	ImgTypeAST=0;	//defines the type of image: 0 - byte, 1  short, 2 - ushort, 3 - int, 4 - long, 5 - ulong, 6 - float, 7- double

	TOF =0;			//time of flight, in seconds
	TimeBin=0;		//width of the time bin, in seconds
	Counts=0;		//counts in the image
	Triggers=0;		//number of triggers received by that image
	ShtrN=0;		//shutter number at which the image was acquired, i.e. TOF belongs to shtrN
}

/*--------------------------------------------------------------------*/
FitsFile::~FitsFile(void)
{
	#ifdef DEBUG_LOG_FILE
	if( _fLog.is_open() ) _fLog.close(); //close the log file
	#endif

	if(FileName != NULL) delete [] FileName;
	if(Img != NULL)		 delete [] Img;
	FileName=NULL;		
	iX=0;				
	iY=0;				
}

/*--------------------------------------------------------------------*/
// ImgType_par defines the type of image: 0 - byte, 1  short, 2 - ushort, 3 - int, 4 - long, 5 - ulong, 6 - float, 7- double
int FitsFile::CreateEmptyImage(int ImgType_par, int iX_par, int iY_par) {
 int pixBytes;

	iX = iX_par, iY = iY_par; ImgTypeAST = ImgType_par;

    pixBytes = getImgBytesPerPixel(ImgType_par);

	Img = new char [iX*iY*pixBytes];
	ZeroImage(); //sets entire image to 0 value

	return 0;
}

/*--------------------------------------------------------------------*/
int FitsFile::CopyImg(FitsFile *ImgToCopy) {
 double dTmp;

	//initialize image if it was not initialized yet
	if(Img == NULL) CreateEmptyImage(ImgToCopy->ImgTypeAST, ImgToCopy->iX, ImgToCopy->iY);

	if(iX != ImgToCopy->iX || iY != ImgToCopy->iY) {LOGCALL("Err copying image: wrong dimensions"); return FITS_COPYIMG_ERROR;}

	TOF = ImgToCopy->TOF;
	TimeBin = ImgToCopy->TimeBin;
	Counts = ImgToCopy->Counts;
	Triggers = ImgToCopy->Triggers;
	ShtrN = ImgToCopy->ShtrN;
	FileIndex = ImgToCopy->FileIndex;

	if(ImgTypeAST == ImgToCopy->ImgTypeAST)  {
		//can just copy entire array into new image
		int bytesPerPixel = getImgBytesPerPixel(ImgTypeAST);
		memcpy_s(Img, iX*iY*bytesPerPixel, ImgToCopy->Img, ImgToCopy->iX*ImgToCopy->iY*bytesPerPixel);
	}
	else {
		for (long i=0; i<iY*iX; i++) {
		    switch (ImgToCopy->ImgTypeAST) {
				case 0:  dTmp = ((char *)ImgToCopy->Img)[i];   break;
				case 1:  dTmp = ((short *)ImgToCopy->Img)[i];  break;
				case 2:  dTmp = ((short *)ImgToCopy->Img)[i];  break;
				case 3:  dTmp = ((int *)ImgToCopy->Img)[i];    break;
				case 4:  dTmp = ((long *)ImgToCopy->Img)[i];   break;
				case 5:  dTmp = ((long *)ImgToCopy->Img)[i];   break;
				case 6:  dTmp = ((float *)ImgToCopy->Img)[i];  break;
				case 7:  dTmp = ((double *)ImgToCopy->Img)[i]; break;
				default: LOGCALL("Err copying image: wrong Img type"); return FITS_COPYIMG_ERROR; break;
			}
			switch (ImgTypeAST) {
				case 0:  ((char *)Img)[i] = (char)dTmp; break;
				case 1:  ((short *)Img)[i] = (short)dTmp; break;
				case 2:  ((short *)Img)[i] = (short)dTmp; break;
				case 3:  ((int *)Img)[i] = (int)dTmp; break;
				case 4:  ((long *)Img)[i] = (long)dTmp; break;
				case 5:  ((long *)Img)[i] = (long)dTmp; break;
				case 6:  ((float *)Img)[i] = (float)dTmp; break;
				case 7:  ((double *)Img)[i] = (double)dTmp; break;
				default: LOGCALL("Err copying image: wrong Img type"); return FITS_COPYIMG_ERROR; break;
			}
		}
	}

	return 0;
}

/*--------------------------------------------------------------------*/
int FitsFile::AddImg(FitsFile *ImgToAdd) {
 double dTmp;

	if(iX != ImgToAdd->iX || iY != ImgToAdd->iY || Img == NULL) {LOGCALL("Err adding image: wrong dimensions or empty image"); return FITS_ADDIMG_ERROR;}

	Counts = 0;

		for (long i=0; i<iY*iX; i++) {
		    switch (ImgToAdd->ImgTypeAST) {
				case 0:  dTmp = ((char *)ImgToAdd->Img)[i];   break;
				case 1:  dTmp = ((short *)ImgToAdd->Img)[i];  break;
				case 2:  dTmp = ((short *)ImgToAdd->Img)[i];  break;
				case 3:  dTmp = ((int *)ImgToAdd->Img)[i];    break;
				case 4:  dTmp = ((long *)ImgToAdd->Img)[i];   break;
				case 5:  dTmp = ((long *)ImgToAdd->Img)[i];   break;
				case 6:  dTmp = ((float *)ImgToAdd->Img)[i];  break;
				case 7:  dTmp = ((double *)ImgToAdd->Img)[i]; break;
				default: LOGCALL("Err adding image: wrong Img type"); return FITS_ADDIMG_ERROR; break;
			}
			switch (ImgTypeAST) {
				case 0:  ((char *)Img)[i] += (char)dTmp; Counts += ((char *)Img)[i]; break;
				case 1:  ((short *)Img)[i] += (short)dTmp; Counts += ((short *)Img)[i]; break;
				case 2:  ((short *)Img)[i] += (short)dTmp; Counts += ((short *)Img)[i]; break;
				case 3:  ((int *)Img)[i] += (int)dTmp; Counts += ((int *)Img)[i]; break;
				case 4:  ((long *)Img)[i] += (long)dTmp; Counts += ((long *)Img)[i]; break;
				case 5:  ((long *)Img)[i] += (long)dTmp; Counts += ((long *)Img)[i]; break;
				case 6:  ((float *)Img)[i] += (float)dTmp; Counts += (unsigned long)((float *)Img)[i]; break;
				case 7:  ((double *)Img)[i] += (double)dTmp; Counts += (unsigned long)((double *)Img)[i]; break;
				default: LOGCALL("Err adding image: wrong Img type"); return FITS_ADDIMG_ERROR; break;
			}
		}

	return 0;
}

/*--------------------------------------------------------------------*/
int FitsFile::DivideByImage(FitsFile *ImgDivideBy) {
 double dTmp, dCounts;

	if(iX != ImgDivideBy->iX || iY != ImgDivideBy->iY || Img == NULL) {LOGCALL("Err dividing image: wrong dimensions or empty image"); return FITS_DIVIMG_ERROR;}

	dCounts = 0;

	for (long i=0; i<iY*iX; i++) {
	    switch (ImgDivideBy->ImgTypeAST) {
			case 0:  dTmp = ((char *)ImgDivideBy->Img)[i];   break;
			case 1:  dTmp = ((short *)ImgDivideBy->Img)[i];  break;
			case 2:  dTmp = ((short *)ImgDivideBy->Img)[i];  break;
			case 3:  dTmp = ((int *)ImgDivideBy->Img)[i];    break;
			case 4:  dTmp = ((long *)ImgDivideBy->Img)[i];   break;
			case 5:  dTmp = ((long *)ImgDivideBy->Img)[i];   break;
			case 6:  dTmp = ((float *)ImgDivideBy->Img)[i];  break;
			case 7:  dTmp = ((double *)ImgDivideBy->Img)[i]; break;
			default: LOGCALL("Err adding image: wrong Img type"); return FITS_ADDIMG_ERROR; break;
		}
		if(dTmp) { //do not divide by 0
			switch (ImgTypeAST) {
				case 0:  ((char *)Img)[i] = (char)(((char *)Img)[i]/dTmp); dCounts += ((char *)Img)[i]; break;
				case 1:  ((short *)Img)[i] = (short)(((short *)Img)[i]/dTmp); dCounts += ((short *)Img)[i]; break;
				case 2:  ((short *)Img)[i] = (short)(((short *)Img)[i]/dTmp); dCounts += ((short *)Img)[i]; break;
				case 3:  ((int *)Img)[i] = (int)(((int *)Img)[i]/dTmp); dCounts += ((int *)Img)[i]; break;
				case 4:  ((long *)Img)[i] = (long)(((long *)Img)[i]/dTmp); dCounts += ((long *)Img)[i]; break;
				case 5:  ((long *)Img)[i] = (long)(((long *)Img)[i]/dTmp); dCounts += ((long *)Img)[i]; break;
				case 6:  ((float *)Img)[i] = (float)(((float *)Img)[i]/dTmp); dCounts += ((float *)Img)[i]; break;
				case 7:  ((double *)Img)[i] = (double)(((double *)Img)[i]/dTmp); dCounts += ((double *)Img)[i]; break;
				default: LOGCALL("Err adding image: wrong Img type"); return FITS_DIVIMG_ERROR; break;
			}
		}
	}
	Counts = (unsigned long)dCounts;
	return 0;
}

/*--------------------------------------------------------------------*/
//int ReadImage_FITS(char *filename, short *Img, unsigned long *Ntiggers, unsigned long *Ncounts, float *TimeBin, float *TOF)
/*--------------------------------------------------------------------*/
int FitsFile::ReadFitsFile(char *FileName_par)
{
 fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
 int status, anynull, iImgDim, imgType_local, pixBytes;
 long naxes[2], StartPixel, npixels;
 float nullval;
 char Comment[256], sTmp[32];
 int i = 0, i1 = 0, i2 = 0, j = 0;

// unsigned long Ntiggers, Ncounts;
// float TimeBin, TOF;
 int dataType;

    status = 0;

    fits_open_file(&fptr, FileName_par, READONLY, &status);
	if (status) return FITS_READ_ERROR;

	//Read the image dimension
	fits_get_img_dim(fptr, &iImgDim, &status);
	if (status || iImgDim != 2) return FITS_READ_ERROR;

	//Read the image dimension
	fits_get_img_size(fptr, iImgDim, naxes, &status); 
	if (status) return FITS_READ_ERROR;
	iX = naxes[0]; iY = naxes[1];

	//Read the image type
	fits_get_img_type(fptr, &imgType_local, &status);
	if (status) return FITS_READ_ERROR;

    /* Create the primary array image of a given type */
    switch (imgType_local) {
		case BYTE_IMG:	 ImgTypeAST = 0;  dataType=TBYTE;   pixBytes=sizeof(char); break;
		case SHORT_IMG:	 ImgTypeAST = 1;  dataType=TSHORT;  pixBytes=sizeof(short); break;
		case LONG_IMG:	 ImgTypeAST = 4;  dataType=TLONG;   pixBytes=sizeof(long); break;
		case FLOAT_IMG:	 ImgTypeAST = 6;  dataType=TFLOAT;  pixBytes=sizeof(float); break;
		case DOUBLE_IMG: ImgTypeAST = 7;  dataType=TDOUBLE; pixBytes=sizeof(double); break;
		default: LOGCALL("Err fits img type"); return FITS_READ_ERROR;  break;
	}

    /* read the keywords */
	if ( fits_read_key(fptr, TLONG, "N_COUNTS",  &Counts,  Comment, &status) ) {Counts=0; status=0;LOGCALL1(FileName_par," No keyword N_COUNTS");}
	if ( fits_read_key(fptr, TLONG, "N_TRIGS",   &Triggers, Comment, &status) ) {Triggers=0;status=0;LOGCALL1(FileName_par," No keyword N_TRIGS");}
	if ( fits_read_key(fptr, TFLOAT, "TIMEBIN",  &TimeBin,  Comment, &status) ) {TimeBin=0; status=0;LOGCALL1(FileName_par," No keyword TIMEBIN");}
	if ( fits_read_key(fptr, TFLOAT, "TOF",      &TOF,      Comment, &status) ) {TOF=0;	    status=0;LOGCALL1(FileName_par," No keyword TOF");}
	if ( fits_read_key(fptr, TFLOAT, "SHUTTER_N",&ShtrN,    Comment, &status) ) {ShtrN=0;   status=0;/*LOGCALL1(FileName_par," No keyword SHUTTER_N");*/}

//	if ( fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) ) {fits_report_error(stderr, status );  return NO_DATA;}

    npixels  = naxes[0] * naxes[1];         /* number of pixels in the image */
	Img = new char [pixBytes*npixels];
	StartPixel = 1; //start reading from the first pixel in image
	fits_read_img(fptr, dataType, StartPixel, npixels, &nullval, Img, &anynull, &status);
	if (status) return FITS_READ_ERROR;

	if ( fits_close_file(fptr, &status) ) {return FITS_READ_ERROR;}

	size_t len = strnlen(FileName_par, 2048);
	FileName = new char [len+1];
	sprintf_s(FileName,len+1,"%s",FileName_par);

	//find the file name index here
	for (int i=len; i>=0; i--) {
		if( FileName[i] == '.') i2=i-1;
		else if ( FileName[i] == '_') { i1=i+1; break;}
	}
	for (i=i1, j=0; i<=i2; i++) { sTmp[j++] = FileName[i];}
	sTmp[j]='\0';
	sscanf_s(sTmp,"%d",&FileIndex);

    return status;
}

//********************************************************************************************
//		Save image into fits file
//		ImageType defines the type of image: 0 - byte, 1  short, 2 - ushort, 3 - int, 4 - long, 5 - ulong, 6 - float, 7- double
//		TBYTE, TSHORT, TUSHORT, TINT, TLONG, TULONG, TFLOAT, TDOUBLE.
//		---------------------------------------------------------------			
//		BYTE_IMG   =   8   ( 8-bit byte pixels, 0 - 255)
//		SHORT_IMG  =  16   (16 bit integer pixels)
//		LONG_IMG   =  32   (32-bit integer pixels)
//		FLOAT_IMG  = -32   (32-bit floating point pixels)
//		DOUBLE_IMG = -64   (64-bit floating point pixels)
//********************************************************************************************
int FitsFile::SaveImage_to_FITS(char *FileName_par) 
{
 fitsfile *fptr;       /* pointer to the FITS file; defined in fitsio.h */
 int status;
 long  fpixel = 1, naxis = 2, nelements;
 long naxes[2];   /* image is 2 dimensional*/
 int fileType, dataType;
 double dTMP;
 unsigned short usTMP;
 unsigned long ulTMP;
 long lTMP;
 char *FileNameLocal;

    naxes[0] = iX; naxes[1] = iY; /* image is iX pixels wide by iY rows */

    status = 0;         /* initialize status before calling fitsio routines */

	//save Filename into local variable
	size_t len = strnlen(FileName_par, 2048);
	FileName = new char [len+1];
	sprintf_s(FileName,len+1,"%s",FileName_par);

	//need that local name for the overwrite files, which is done by adding ! in front of the file name. 
	FileNameLocal = new char [len+2];
	sprintf_s(FileNameLocal,len+2,"!%s",FileName_par);

    fits_create_file(&fptr, FileNameLocal, &status);   /* create new file, use that name with ! added in front to overwrite the file*/
	delete [] FileNameLocal;
	if (status) return FITS_WRITE_ERROR;
    //fits_report_error(stderr, status);  /* print out any error messages */

    /* define the variable for saving the file of image type and data type*/
    switch (ImgTypeAST) {
		case 0:  fileType = BYTE_IMG;   dataType=TBYTE;   break;
		case 1:  fileType = SHORT_IMG;  dataType=TSHORT;  break;
		case 2:  fileType = SHORT_IMG;  dataType=TUSHORT; break;
		case 3:  fileType = SHORT_IMG;  dataType=TINT;    break;
		case 4:  fileType = LONG_IMG;   dataType=TLONG;   break;
		case 5:  fileType = LONG_IMG;   dataType=TULONG;  break;
		case 6:  fileType = FLOAT_IMG;  dataType=TFLOAT;  break;
		case 7:  fileType = DOUBLE_IMG; dataType=TDOUBLE; break;
		default: fileType = FLOAT_IMG;  dataType=TFLOAT;  break;
	}

    /* Create the primary array image of type set by variable ImgType*/
    fits_create_img(fptr, fileType, naxis, naxes, &status);
	if (status) goto EXIT_HERE;
    //fits_report_error(stderr, status);  /* print out any error messages */

	dTMP= TOF;		fits_update_key(fptr, TDOUBLE, "TOF",	   &dTMP, "Ttime of flight from the external trigger", &status);
	dTMP= TimeBin;	fits_update_key(fptr, TDOUBLE, "TIMEBIN",  &dTMP, "Time width of this image", &status);
	ulTMP = Counts;	fits_update_key(fptr, TULONG,  "N_COUNTS", &ulTMP, "Total counts in this image", &status);
	lTMP = Triggers;fits_update_key(fptr, TLONG,   "N_TRIGS",  &lTMP, "Number of triggers acquired", &status);
	usTMP = ShtrN;	fits_update_key(fptr, TUSHORT, "SHUTTER_N", &usTMP, "The shutter number at which this image was acquired, starts from 1", &status);

    nelements = naxes[0] * naxes[1];          /* number of pixels to write */

    /* Write the array of integers to the image */
    fits_write_img(fptr, dataType, fpixel, nelements, (void *)Img, &status);
	if (status) goto EXIT_HERE;
	//fits_report_error(stderr, status);  /* print out any error messages */

EXIT_HERE:
    fits_close_file(fptr, &status);            /* close the file */
	if (status) return FITS_WRITE_ERROR;
	//fits_report_error(stderr, status);  /* print out any error messages */

    return( status );
}

/*--------------------------------------------------------------------*/
int FitsFile::getImgBytesPerPixel(int ImgType_par)
{
 int pixBytes;

    /* define the variable for saving the file of image type and data type*/
	//	ImageType defines the type of image: 0 - byte, 1  short, 2 - ushort, 3 - int, 4 - long, 5 - ulong, 6 - float, 7- double
	//	TBYTE, TSHORT, TUSHORT, TINT, TLONG, TULONG, TFLOAT, TDOUBLE.
    switch (ImgType_par) {
		case 0:  pixBytes=sizeof(char);   break;
		case 1:  pixBytes=sizeof(short);  break;
		case 2:  pixBytes=sizeof(short);  break;
		case 3:  pixBytes=sizeof(int);    break;
		case 4:  pixBytes=sizeof(long);   break;
		case 5:  pixBytes=sizeof(long);   break;
		case 6:  pixBytes=sizeof(float);  break;
		case 7:  pixBytes=sizeof(double); break;
		default: pixBytes=sizeof(short);  break;
	}

	return pixBytes;
}

/*--------------------------------------------------------------------*/
// Set all elements of the image to 0
/*--------------------------------------------------------------------*/
int FitsFile::ZeroImage(void) {
 int pixBytes;

	pixBytes = getImgBytesPerPixel(ImgTypeAST);

	memset(Img, 0, iX*iY*pixBytes);
	Counts = 0;

	return 0;
}

/*--------------------------------------------------------------------*/
// Divide image by a value
/*--------------------------------------------------------------------*/
int FitsFile::DivideByConst(double Divider) {
 double dCounts;
 unsigned long nPixs;

	if(iX <= 0 || iY <= 0 || Img == NULL) {LOGCALL("Err dividing image: wrong dimensions or empty image"); return FITS_DIVIMG_ERROR;}
	//do not divide by 0
	if(!Divider)  {LOGCALL("Err dividing image: divider value is 0 "); return FITS_DIVIMG_ERROR;}

	dCounts = 0;
	nPixs = (unsigned long)iY*iX;
	for (unsigned long i=0; i<nPixs; i++) {			
		switch (ImgTypeAST) {
				case 0:  ((char *)Img)[i] = (char)(((char *)Img)[i]/Divider); dCounts += ((char *)Img)[i]; break;
				case 1:  ((short *)Img)[i] = (short)(((short *)Img)[i]/Divider); dCounts += ((short *)Img)[i]; break;
				case 2:  ((short *)Img)[i] = (short)(((short *)Img)[i]/Divider); dCounts += ((short *)Img)[i]; break;
				case 3:  ((int *)Img)[i] = (int)(((int *)Img)[i]/Divider); dCounts += ((int *)Img)[i]; break;
				case 4:  ((long *)Img)[i] = (long)(((long *)Img)[i]/Divider); dCounts += ((long *)Img)[i]; break;
				case 5:  ((long *)Img)[i] = (long)(((long *)Img)[i]/Divider); dCounts += ((long *)Img)[i]; break;
				case 6:  ((float *)Img)[i] = (float)(((float *)Img)[i]/Divider); dCounts += ((float *)Img)[i]; break;
				case 7:  ((double *)Img)[i] = (double)(((double *)Img)[i]/Divider); dCounts += ((double *)Img)[i]; break;
				default: LOGCALL("Err adding image: wrong Img type"); return FITS_DIVIMG_ERROR; break;
		}
	}
	Counts = (unsigned long)dCounts;
	return 0;
}


/*--------------------------------------------------------------------*/
// get number of counts in the given area of the image
/*--------------------------------------------------------------------*/
double FitsFile::getCountsInArea(unsigned short X0, unsigned short X1, unsigned short Y0, unsigned short Y1) {
 double dCounts=0;

	if( X1 > iX || Y1 > iY || X0<0 || Y0<0) {LOGCALL("Cannot calculate counts in the imageL: wrong area defined in parameters"); return FITS_BUILD_SPECTRA_ERR;} 

	for(int k=Y0; k <= Y1; k++) {
		for(int m=X0; m <= X1; m++) {
			switch (ImgTypeAST) {
				case 0:  dCounts += ((char *)Img)  [k*iX + m]; break;
				case 1:  dCounts += ((short *)Img) [k*iX + m]; break;
				case 2:  dCounts += ((short *)Img) [k*iX + m]; break;
				case 3:  dCounts += ((int *)Img)   [k*iX + m]; break;
				case 4:  dCounts += ((long *)Img)  [k*iX + m]; break;
				case 5:  dCounts += ((long *)Img)  [k*iX + m]; break;
				case 6:  dCounts += ((float *)Img) [k*iX + m]; break;
				case 7:  dCounts += ((double *)Img)[k*iX + m]; break;
				default: LOGCALL("Err counts in image: wrong Img type"); return FITS_BUILD_SPECTRA_ERR; break;
			}
		}
	}
	return dCounts;
}
