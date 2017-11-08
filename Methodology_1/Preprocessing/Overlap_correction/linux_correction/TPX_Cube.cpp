#include <iostream>
#include <fstream>	
#include <string.h>
#include <math.h>
#include <libgen.h>
#include "TPX_Cube.h"
#include "../common/build_file_list.h"
#include "../common/error_handling.h"
#include "../common/output_directory.h"
#include "../common/utils.h"

using namespace std;

#define sprintf_s(a, b, c, ...) snprintf(a, b, c, __VA_ARGS__)
#define memcpy_s(a, b, c, d) memcpy(a, c, d)
#define sscanf_s sscanf
#define strcpy_s(a, b, c) strcpy(a, c)


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
TPX_Cube::TPX_Cube(void)
{
	FileNameBase=NULL;	
	Dir_Name=NULL;
	nSlices =0;
	nShtrs =0;
	TOF_Arr=NULL;
	Shtr_Arr= NULL;
	TimeBin_Arr = NULL;
	Wvl_Arr=NULL;
	E_Arr = NULL;
	L=0; Tshift=0;
	Slices=NULL;

	Spectra_Arr=NULL;
	SpectrIndLeft=0; SpectrIndRight=0;
}

/*********************************************************************************************/
//	Create TPX cube from the data in a given directory. Filname - is one file from the data set with full path to it
/*--------------------------------------------------------------------*/
TPX_Cube::TPX_Cube(char *FileName_par)
{
	FileNameBase=NULL;	
	Dir_Name=NULL;
	nSlices =0;
	nShtrs =0;
	TOF_Arr=NULL;
	Shtr_Arr= NULL;
	TimeBin_Arr = NULL;
	Wvl_Arr=NULL;
	E_Arr = NULL;
	L=0; Tshift=0;
	Slices=NULL;

	Spectra_Arr=NULL;
	SpectrIndLeft=0; SpectrIndRight=0;

	LoadFromFiles(FileName_par);
}

/*--------------------------------------------------------------------*/
TPX_Cube::~TPX_Cube(void)
{
	if(Slices != NULL) delete [] Slices;
	if(TOF_Arr != NULL) delete [] TOF_Arr;
	if(Spectra_Arr != NULL) delete [] Spectra_Arr;
	if(Shtr_Arr != NULL) delete [] Shtr_Arr;
	if(TimeBin_Arr != NULL) delete [] TimeBin_Arr;
	if(Wvl_Arr != NULL) delete [] Wvl_Arr;
	if(E_Arr != NULL) delete [] E_Arr;
	if(FileNameBase != NULL) delete [] FileNameBase;
	if(Dir_Name != NULL) delete [] Dir_Name;
	
	FileNameBase=NULL;
	Dir_Name=NULL;
	nSlices =0;
	nShtrs =0;
	TOF_Arr=NULL;
	Shtr_Arr= NULL;
	TimeBin_Arr = NULL;
	Wvl_Arr=NULL;
	E_Arr = NULL;
	L=0; Tshift=0;
	Slices=NULL;

	Spectra_Arr=NULL;
	SpectrIndLeft=0; SpectrIndRight=0;
}

static char* dirname_basename(char* input, bool get_dirname)
{
	char* copy = (char*)calloc(1, strlen(input) + 1);
	strcpy(copy, input);

	char* modified = get_dirname ? dirname(copy) : basename(copy);
	char* tmp = new char[strlen(modified) + 1];
	strcpy(tmp, modified);
	free(copy);
	return tmp;
}

/*static char* dirname_basename(char* input, bool get_dirname)
{
    int i=0, index = -1, len = strlen(input);
    for (i=0;i<len;i++)
    {
        if (input[i] == '/')
        {
            index = i;
        }
    }
    
    char* copy = (char*)calloc(1, len + 1);
    
    if (get_dirname)
    {
        strncpy(copy, input, index);
    }
    else
    {
        index++;
        strcpy(copy, &input[index]);
    }
    
    return copy;
}*/

/*********************************************************************************************/
//	Load TPX Cube from files and set all the parameters of the cube from the directory data
/*-------------------------------------------------------------------------------------------*/
int TPX_Cube::LoadFromFiles(char *FileName_par)
{
 int i,j, ret, iTmp;
 double dTmp1, dTmp2;
 char *FileNameToRead;
 ifstream infile;

	//form the file name base which will be used to read the TPX Cube
	size_t len = strnlen(FileName_par, 2048); if(len<=11) {LOGCALL1("NO file name passed to constructor of TPX Cube",len); return -1;}
	FileNameBase = new char [len+1];
	sprintf_s(FileNameBase,len+1,"%s",FileName_par);
	FileNameBase[len-11]='\0'; //cut off _xxxxx.fits from the file name

	FileNameToRead = new char [len+64]; //buffer for the file name which will be read

	//find out how many files are in that directory - taht will be the max value for the numerb of files to be read in
	char* directoryName = dirname_basename(FileName_par, true);

	char* filename = dirname_basename(FileName_par, false);

    char fileSearchPattern[128] = {0};
	strncpy(fileSearchPattern, filename, strlen(filename) - 11);
	strcat(fileSearchPattern, "*.fits\0");

	FILE_LIST flist;
	ret = build_file_list(&flist, directoryName, fileSearchPattern);
	if (ret != ERH_NO_ERROR)
		return ERROR("couldn't build file list", ret);

	Dir_Name = new char [len+10];
	sprintf_s(Dir_Name,len+10,"%s", directoryName);

	int FilesN = flist.count;
	FilesN -= 1; //minus one as the file SummedImg should be excluded

	//allocate memory for Images in Cube
	Slices = new FitsFile [FilesN]; 
	TOF_Arr = new double [FilesN]; 
	Spectra_Arr = new double [FilesN]; 
	Shtr_Arr = new unsigned short [FilesN]; 
	TimeBin_Arr = new double [FilesN]; 
	nSlices = 0;

	//read Images from fits files
	for(i=0; i<FilesN; i++) {
		sprintf_s(FileNameToRead, len+1, "%s_%05d.fits",FileNameBase,i);
		ret = Slices[i].ReadFitsFile(FileNameToRead);
		if(ret==0) {
			nSlices++;		//number of spectral elements
			TOF_Arr[i]=Slices[i].TOF;		
			Spectra_Arr[i] = Slices[i].Counts;
			TimeBin_Arr[i] = Slices[i].TimeBin;
		}
	}

	//read shutter time values
	//that section was for using file ShuterTimes, but it does not work as the error is cummulitive. 
	//Need to correct it TimeBin_Arr or by Sptra.txt which is the same array
	sprintf_s(FileNameToRead, len+64, "%s_ShutterTimes.txt",FileNameBase);
	infile.open(FileNameToRead); 
	if(infile.fail()) LOGCALL1("Shutter times file was not found",FileNameToRead);
	else {
		i=0;
		while (!infile.eof() && i<MAX_SHUTTERS_AVAILABLE)   {
			infile >> iTmp >> dTmp1 >> dTmp2;
			ShtrStart[i] = dTmp1; 
			if (i>0 && dTmp2>0) ShtrStart[i] += ShtrEnd[i-1];
			ShtrEnd[i] = ShtrStart[i] + dTmp2;
			if( ShtrEnd[i] || ShtrStart[i]) { //do correction only when they are real values, not all the 0,0 at the end of each trigger
			//get correction here to the nearest value of TimeBin_Arr, the values in ShutterTimes.txt are rounded and lead to cummulutive error
				ShtrStart[i] = GetNearestArrValue(ShtrStart[i], TOF_Arr,  FilesN);
				ShtrEnd[i]   = GetNearestArrValue(ShtrEnd[i], TOF_Arr,  FilesN);
			}
			i++;
		}
	}
	if( infile.is_open() ) infile.close();  //close the file

	//find which shutter each slice belongs to
	for(i=0; i<FilesN; i++) {
		for (j=0; j<MAX_SHUTTERS_AVAILABLE; j++) {
			if((TOF_Arr[i] >= ShtrStart[j]-1e-6) && (TOF_Arr[i] <= ShtrEnd[j]+1e-6) ) { 
				Shtr_Arr[i] = j+1; 
				Slices[i].ShtrN=j+1;
				nShtrs = j+1;
				break;
			}
			else if(j==MAX_SHUTTERS_AVAILABLE-1) 
				LOGCALL1("Could not find shutter time for this slice",TOF_Arr[i]);
		}
	}

	//read shutter counts values
	sprintf_s(FileNameToRead, len+64, "%s_ShutterCount.txt",FileNameBase);
	infile.open(FileNameToRead); 
	if(infile.fail()) LOGCALL1("Shutter counts file was not found",FileNameToRead);
	else {
		i=0;
		while (!infile.eof() && i<MAX_SHUTTERS_AVAILABLE)   infile >> j >> ShtrCounts[i++];
	}
	if( infile.is_open() ) infile.close();  //close the file

	SpectrIndLeft = 0;	SpectrIndRight = nSlices-1; //current indexes of the spectra area - Full spectra
	SpectrX0=0; SpectrX1=Slices[0].iX; SpectrY0=0; SpectrY1=Slices[0].iY;		 //current position of the spectra area - full image

	LOGCALL2("Read ", nSlices," images in the TPX cube");

	delete [] FileNameToRead;

	return 0;
}

/*********************************************************************************************/
//	Save TPX Cube into FITS files
/*-------------------------------------------------------------------------------------------*/
int TPX_Cube::SaveCubeToFitsFiles(char *FileNameBase_par)
{
	size_t len = strnlen(FileNameBase_par, 2048);
	char *FileNameLocal = new char [len+64];

	for(int i=0; i<nSlices; i++) {
		sprintf_s(FileNameLocal,len+64,"%s_%05d.fits",FileNameBase_par,Slices[i].FileIndex);
		Slices[i].SaveImage_to_FITS(FileNameLocal);
		//assign new name to that slice
		if(Slices[i].FileName != NULL) delete [] Slices[i].FileName;
		size_t len = strnlen(FileNameLocal, 2048);
		Slices[i].FileName = new char [len+1];
		strcpy_s(Slices[i].FileName, len+1, FileNameLocal);
	}
	return 0;
}

/*********************************************************************************************/
//	create copy of existing Cube, but the image file type can be changed - need that for the floating point of normalized of overlap corrected cube
//	0 - byte, 1  short, 2 - ushort, 3 - int, 4 - long, 5 - ulong, 6 - float, 7- double
/*--------------------------------------------------------------------------------------------*/
int	TPX_Cube::CreateCopy(int ImgType_par, TPX_Cube *SrcCube) {
 FitsFile *SrcSlice;

	nSlices = SrcCube->nSlices;		//number of spectral elements
	nShtrs = SrcCube->nShtrs;		//number of shutters used for acquisition
	for (int i=0; i<MAX_SHUTTERS_AVAILABLE;i++) {
		ShtrCounts[i]=SrcCube->ShtrCounts[i];//array of count values in each of 255 chutters available
		ShtrStart[i]=SrcCube->ShtrStart[i];  //array of start times for shutters
		ShtrEnd[i]=SrcCube->ShtrEnd[i];		 //array of end times for shutters
	}

	Slices = new FitsFile [nSlices]; 
	TOF_Arr = new double [nSlices]; 
	Spectra_Arr = new double [nSlices]; 
	Shtr_Arr = new unsigned short [nSlices]; 
	TimeBin_Arr = new double [nSlices]; 

	size_t len = strnlen(SrcCube->FileNameBase, 2048);
	FileNameBase = new char [len+1];
	strcpy_s(FileNameBase, len+1, SrcCube->FileNameBase);
	Dir_Name	 = new char [len+10];
	strcpy_s(Dir_Name, len+10, SrcCube->Dir_Name);

	for(int i=0; i<nSlices; i++) {
		SrcSlice = &(SrcCube->Slices[i]);
		Slices[i].CreateEmptyImage(ImgType_par, SrcSlice->iX, SrcSlice->iY); 
		Slices[i].CopyImg(SrcSlice); 
		TOF_Arr[i] = SrcSlice->TOF;	//array of TOF values
		Shtr_Arr[i] = SrcSlice->ShtrN;	//array of shutter number for each slice 
		TimeBin_Arr[i] = SrcSlice->TimeBin;	//array of time bin values for each slice
		size_t len = strnlen(SrcSlice->FileName, 2048);
		Slices[i].FileName = new char [len+1];
		strcpy_s(Slices[i].FileName, len+1, SrcSlice->FileName);
		Slices[i].FileIndex = SrcSlice->FileIndex;
	}
	return 0;
}

/*********************************************************************************************/
//	create Images with overlap correction and save them into new Cube file set, same name indexes
//	the new files will be floating point as counts will be non-integers now
/*-------------------------------------------------------------------------------------------*/
int	TPX_Cube::SaveOverlapCorrectedCube(char *FileNameBase_par, bool ROTAX){
 FitsFile *CntsIntgr, *CorrectedFrame;
 int CurrentShtrN=-1;
 double dCounts, dShtrNRatio;
 float *CorrImg; 
 unsigned long *CntsImg;
 unsigned long nPixs, j;
 double nTrigs;
 char sTmp[1024];
 int ret = 0;

	if(nSlices <= 0) return CUBE_N_SLICES_ERR;

	CntsIntgr = new  FitsFile;
	CorrectedFrame = new  FitsFile;
	
	CntsIntgr->CreateEmptyImage(5 /* ulong type image */, Slices[0].iX, Slices[0].iY);
	CorrectedFrame->CreateEmptyImage(6 /* float type image */, Slices[0].iX, Slices[0].iY);
	CorrImg = (float *)(CorrectedFrame->Img);
	CntsImg = (unsigned long *)(CntsIntgr->Img);

	//create subdirectory for the files to be saved
	char* directoryName = dirname_basename(FileNameBase_par, true);
	directoryName = path_join(directoryName, "Corrected");
	ret = mkdir_output_folder(directoryName);
	if (ret != ERH_NO_ERROR)
		return ERROR("couldn't make output folder", ret);

	char* FileNameRoot = dirname_basename(FileNameBase_par, false);

	for(int i=0; i<nSlices; i++) {
		//restart the occupancy matrix if needed - start of new shutter condition
		if(CurrentShtrN < Shtr_Arr[i]) { 
			CntsIntgr->ZeroImage(); 
			CurrentShtrN = Shtr_Arr[i];
			if(CurrentShtrN > nShtrs+1) {LOGCALL2("Too many shutter interrupts found in overlap correction, more than used in the measurements!",CurrentShtrN, nShtrs); ret=CUBE_OVERLAP_CORR_ERR; goto ExitHere;}
			nTrigs = ShtrCounts[CurrentShtrN-1]; //shutter number starts from 1, that is why -1
			dShtrNRatio = (double)nTrigs / ShtrCounts[0]; //normalization factor for the ratio between the shutters as they may hae different number of triggers
		}

		//create floating point copy of the slice to be corrected
		CorrectedFrame->CopyImg(&(Slices[i]));

		//calculate the corrected matrix now
		nPixs = (unsigned long)Slices[i].iY*Slices[i].iX;
		if(CntsIntgr->Counts && nTrigs>0 ) { //there are non-zero cells in the frame and nTriggers is larger than 0
			dCounts = 0;
			for (j=0; j<nPixs; j++) {			
				if(CntsImg[j] > nTrigs) { // that pixel has more counts than triggres which can happen when the value exceeds the maximum short integer number of ~65000 or comes as negative value in IMageJ
					LOGCALL2("Too many counts in occupancy matrix, more than triggers. Can only have one count per trigger!!",CntsImg[j],nTrigs); 
					CorrImg[j] = 0; //make that value to 0, that pixel is not giving us any reasonable data, too noisy!
					//ret=CUBE_OVERLAP_CORR_ERR; 
					//goto ExitHere;
				}
				if (ROTAX) CorrImg[j] /= (float)((1.0-(double)CntsImg[j]/(nTrigs*4/5))*dShtrNRatio);   //correction for the ROTAX trigger
				else       CorrImg[j] /= (float)((1.0-(double)CntsImg[j]/(nTrigs    ))*dShtrNRatio);   //not ROTAX
				//CorrImg[j] /= (float)((1.0-(double)CntsImg[j]/nTrigs)*dShtrNRatio); 
				dCounts += CorrImg[j];
			}
			CorrectedFrame->Counts = (unsigned long)dCounts;
		}

		//add the current frame to the occupancy matrix
		CntsIntgr->AddImg(&(Slices[i]));

		//save the corrected frame into a fits file
		sprintf_s(sTmp, 1024,"%s/%s_%05d.fits",directoryName,FileNameRoot,Slices[i].FileIndex);
		CorrectedFrame->SaveImage_to_FITS(sTmp);

printf("%s\tCounts: %9lu --> %9lu\n", sTmp, Slices[i].Counts, CorrectedFrame->Counts);
	}

ExitHere:
	delete 	CntsIntgr;
	delete	CorrectedFrame;

	return ret;
}

/*********************************************************************************************/
//	calculate the Wvlgth for a given TOF value, FlightPath needs to be set as well
//	meters and seconds in parameters, returns wavelength in Angtroms
/*--------------------------------------------------------------------*/
double	TPX_Cube::TOF_To_Wvlngth(double TOF, double FlightPath_par, double Tshift_par) {  
	return (TOF + Tshift_par) * 3956 / FlightPath_par;  
}

/*********************************************************************************************/
//	calculate the TOF as it would be measured by the system wiht Tshift delay on trigger and L flight path
//	meters and seconds in L and dT parameters, Wvlngth in angstroms, returns TOF is seconds
/*-------------------------------------------------------------------------------------------*/
double	TPX_Cube::Wvlngth_to_TOF(double Wvlngth, double FlightPath_par, double Tshift_par) {  
	return (Wvlngth / 3956 * FlightPath_par) - Tshift_par;
}

/*********************************************************************************************/
//	calculate the Wvlgth for a given TOF value, FlightPath needs to be set as well
//	meters and seconds in parameters, returns wavelength in Angtroms
/*-------------------------------------------------------------------------------------------*/
double	TPX_Cube::TOF_To_E(double TOF_par, double FlightPath_par, double Tshift_par) {  
 double dTmp;
	dTmp = TOF_To_Wvlngth(TOF_par, FlightPath_par, Tshift_par);
	return (6.626*6.626)/2/1.675/dTmp/dTmp/160;
}

/*********************************************************************************************/
//	calculate the TOF as it would be measured by the system wiht Tshift delay on trigger and L flight path
//	meters and seconds in L and dT parameters, Wvlngth in angstroms, returns TOF is seconds
/*-------------------------------------------------------------------------------------------*/
double	TPX_Cube::E_to_TOF(double E, double FlightPath_par, double Tshift_par) {  
 double dTmp;
	dTmp = sqrt((6.626*6.626)/2/1.675/E/160);
	return dTmp / 3956 * FlightPath_par -Tshift_par;
}

/*********************************************************************************************/
//	sets internal variables for the Flight path and dT shift of the trigger in the experiment
//	then also recalculates the Wvlngth and Energy arrays for the neutrons
//	meters and seconds in L and dT parameters, Wvlngth in angstroms, E in eV
/*-------------------------------------------------------------------------------------------*/
int	TPX_Cube::SetExperimVariables(double FlightPath_par, double Tshift_par) {

	if(FlightPath_par <=0 ) {LOGCALL("The FLight path cannot be 0 or negative"); return CUBE_FLIGHT_PATH_ERR;} 

	L = FlightPath_par;
	Tshift = Tshift_par;

	if(nSlices>0 && TOF_Arr != NULL) {
		if(Wvl_Arr != NULL) delete [] Wvl_Arr;
		if(E_Arr != NULL) delete [] E_Arr;
		Wvl_Arr =  new double [nSlices];
		E_Arr =  new double [nSlices];
		for(int i=0; i<nSlices; i++) {
			Wvl_Arr[i] = TOF_To_Wvlngth(TOF_Arr[i], L, Tshift);
			E_Arr[i] = TOF_To_E(TOF_Arr[i], L, Tshift);
		}
	}
	else {LOGCALL("Set flight path before cube was initialized with real values"); return 0;} 

	return 0;
}


/*********************************************************************************************/
//	Returns spectra for a given area of the image and given spectral range
//	X0, X1, Y0, Y1 are in pixels starting from 0, TOF0 and TOF1 are in Time of flight, which can be not exactly values in the TOF_array
/*-------------------------------------------------------------------------------------------*/
double *TPX_Cube::GetSpectra(unsigned short *iDim, unsigned short X0, unsigned short X1, unsigned short Y0, unsigned short Y1, double TOF0, double TOF1) {
 int iIndLeft=-1, iIndRight=-1;

	if( nSlices <= 0	   || 
		(X1<X0) || (Y1<Y0) || 
		//(X0<0)  || (Y0<0)  || 
		//(X1>Slices[0].iX)  || (Y1>Slices[0].iY) ||
		(TOF1<TOF0)		   || 
		TOF0<TOF_Arr[0]    || TOF1>TOF_Arr[nSlices-1]) 
		{LOGCALL("Cannot build spectra: called with wrong parameters"); return NULL;} 

	//first calculate how many spectral elements we will need to get the spectra between TOF0 and TOF1
	for(int i=0; i<nSlices-1; i++) {
		if(TOF0 >= TOF_Arr[i] && TOF0 < TOF_Arr[i+1]) { iIndLeft=i;	break; }
	}
	if(iIndLeft == -1) {LOGCALL("Cannot find TOF value in the current Cube data: build spectra called with wrong parameters"); return NULL;} 
	for(int i=iIndLeft; i<nSlices-1; i++) {
		if(TOF1 >= TOF_Arr[i] && TOF1 < TOF_Arr[i+1]) { iIndRight=i; break; }
	}
	if(iIndRight == -1) {LOGCALL("Cannot find TOF value in the current Cube data: build spectra called with wrong parameters"); return NULL;} 

	*iDim = iIndRight - iIndLeft + 1;

	for(int i=iIndLeft; i<=iIndRight; i++) Spectra_Arr[i] = Slices[i].getCountsInArea(X0,X1,Y0,Y1);

	SpectrIndLeft = iIndLeft;	SpectrIndRight = iIndRight; //current indexes of the spectra area
	SpectrX0=X0; SpectrX1=X1; SpectrY0=Y0; SpectrY1=Y1;		//current position of the spectra area

	return &(Spectra_Arr[SpectrIndLeft]);
}

/*********************************************************************************************/
//	Saves the last spectra for a given area of the image
/*-------------------------------------------------------------------------------------------*/
int TPX_Cube::SaveLastSpectra(char *FileName){
 ofstream	_fOut;

 if( nSlices <= 0 || SpectrIndLeft<0 || SpectrIndRight<=0) {LOGCALL("Err saving spectral file: No spectra was calculated yet"); return -1;} 

	_fOut.open(FileName);

	if( _fOut.is_open() ) {
		_fOut <<"Filename Base\t"<<FileNameBase<<endl;
		_fOut <<"X0\t"<<SpectrX0<<"\tX1\t"<<SpectrX1<<"\tY0\t"<<SpectrY0<<"\tY1\t"<<SpectrY1<<"\tL\t"<<L<<"\tTshift\t"<<Tshift<<endl;
		_fOut <<"TOF\t"<<"Wavelngth\t"<<"E\t"<<"Counts"<<endl;
			
		for(int i=SpectrIndLeft; i<=SpectrIndRight; i++) 
			_fOut << TOF_Arr[i] << "\t" << TOF_To_Wvlngth(TOF_Arr[i],L,Tshift) << "\t"  << TOF_To_E(TOF_Arr[i],L,Tshift) << "\t" << Spectra_Arr[i] <<endl;

		_fOut.close(); //close the log file
	}
	return 0;
}


/*********************************************************************************************/
//	Finds the nearest value in the array
/*-------------------------------------------------------------------------------------------*/
double TPX_Cube::GetNearestArrValue(double Val, double *Arr, int iDim) {
 double dMin=1e10, dTmp;
 double retValue=0;

	for (int i=0; i< iDim; i++) {
		dTmp = fabs(Val-Arr[i]);
		if(  dTmp < dMin ) {
			dMin = dTmp; 
			retValue = Arr[i];
		}
	}

	if(dMin < 1e10)	return retValue;
	else {
		LOGCALL("Could not find value close to original value in the given array");
		return Val;
	}
}
