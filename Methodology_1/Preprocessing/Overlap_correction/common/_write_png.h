#ifndef WRITE_PNG_H
#define WRITE_PNG_H

#include <png.h>

//----------------------------------------------
//    function prototypes
//----------------------------------------------
int write_png_greyscale16(char* path, int xdim, int ydim, uint16_t* img_data, int nentries, png_text* metadata);

#endif

