#include "utils.h"
#include "error_handling.h"
#include "write_png.h"

int write_png_greyscale16(char* path, int xdim, int ydim, uint16_t* img_data, int nentries, png_text* metadata)
{
	//initialize stuff
	png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	if (png_ptr == NULL)
		return ERROR("[write_png_file] png_create_write_struct failed", ERH_ERROR);

	png_infop info_ptr = png_create_info_struct(png_ptr);
	if (info_ptr == NULL)
		return ERROR("[write_png_file] png_create_info_struct failed", ERH_ERROR);

	//if (setjmp(png_jmpbuf(png_ptr)))
	//	return ERROR("[write_png_file] Error during init_io", ERH_ERROR);

	FILE* fp = fopen(path, "wb");
	if (fp == NULL)
		return ERROR("couldn't open file for writing", ERH_ERROR);

	png_init_io(png_ptr, fp);

	if (nentries > 0)
		png_set_text(png_ptr, info_ptr, metadata, nentries);

	//write header
	png_set_IHDR(	png_ptr,
			info_ptr,
			xdim,
			ydim,
			16,	//bit depth
			PNG_COLOR_TYPE_GRAY,
			PNG_INTERLACE_NONE,
			PNG_COMPRESSION_TYPE_BASE,
			PNG_FILTER_TYPE_BASE);

	png_write_info(png_ptr, info_ptr);

	//big-endian/little-endian correction
	png_set_swap(png_ptr);

	png_bytep* row_pointer = (png_bytep*)malloc(sizeof(png_bytep) * ydim);
	//todo: error checking

	int i = 0;
	for (i=0;i<ydim;i++)
		row_pointer[i] = (uint8_t*)&img_data[i * xdim];

	//write bytes
	png_write_image(png_ptr, row_pointer);

	//end write
	png_write_end(png_ptr, NULL);

	png_destroy_write_struct(&png_ptr, &info_ptr);
	free(row_pointer);
	fclose(fp);
	return ERH_NO_ERROR;
}

