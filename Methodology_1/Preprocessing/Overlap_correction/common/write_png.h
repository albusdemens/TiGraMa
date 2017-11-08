#ifndef WRITE_PNG_H
#define WRITE_PNG_H

typedef struct
{
	uint8_t r;
	uint8_t g;
	uint8_t b;
} RGB;

//----------------------------------------------
//    function prototypes
//----------------------------------------------
int write_png_rgb(char* path, int xdim, int ydim, RGB* rgb_data, int zoom_factor);
int write_png_greyscale16(char* path, int xdim, int ydim, uint16_t* img_data, int zoom_factor);
int write_png_greyscale8(char* path, int xdim, int ydim, uint8_t* img_data, int zoom_factor);
int write_png_binary(char* path, int xdim, int ydim, uint8_t* img_data, int zoom_factor);

#endif

