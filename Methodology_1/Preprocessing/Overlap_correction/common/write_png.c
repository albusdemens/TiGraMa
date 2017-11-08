#include <png.h>
#include "utils.h"
#include "error_handling.h"
#include "write_png.h"

static int write_png_file(char* filepath, uint8_t* buf, size_t bytes_per_row, int width, int height, int bit_depth, int color_type)
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

	FILE* fp = fopen(filepath, "wb");
	if (fp == NULL)
		return ERROR("couldn't open file for writing", ERH_ERROR);

	png_init_io(png_ptr, fp);

	//write header
	png_set_IHDR(	png_ptr,
			info_ptr,
			width,
			height,
			bit_depth,
			color_type,
			PNG_INTERLACE_NONE,
			PNG_COMPRESSION_TYPE_BASE,
			PNG_FILTER_TYPE_BASE);

	png_write_info(png_ptr, info_ptr);

	//big-endian/little-endian correction
	png_set_swap(png_ptr);

	png_bytep* row_pointer = (png_bytep*)malloc(sizeof(png_bytep) * height);
	//todo: error checking

	int i = 0;
	for (i=0;i<height;i++)
		row_pointer[i] = &buf[i * bytes_per_row];

	//write bytes
	png_write_image(png_ptr, row_pointer);

	//end write
	png_write_end(png_ptr, NULL);

	png_destroy_write_struct(&png_ptr, &info_ptr);
	free(row_pointer);
	fclose(fp);
	return ERH_NO_ERROR;
}

static void setb_png(uint8_t* buf, int i)
{
	int wordpos = i / 8;
	int bitpos = 7 - i % 8;

	buf[wordpos] |= 1 << bitpos;
}

static int pack_bits(int xdim, int ydim, int bytes_per_row, uint8_t* expanded, uint8_t* binarized)
{
	int i = 0, j = 0;

	for (j=0;j<ydim;j++)
	{
		for (i=0;i<xdim;i++)
		{
			int v = expanded[j * xdim + i];

			if (v < 0 || v > 1)
				return ERROR("pixel values must be in range [0, 1]", ERH_ERROR);

			if (v == 1)
				setb_png(&binarized[j * bytes_per_row], i);
		}
	}

	return ERH_NO_ERROR;
}

static int zoom_image(int bit_depth, int xdim, int ydim, int zoom_factor, uint8_t* src, uint8_t** p_data)
{
	int i = 0, j = 0;
	int z_xdim = zoom_factor * xdim;
	int z_ydim = zoom_factor * ydim;

	if (bit_depth <= 0 || bit_depth % 8)
		return ERROR("invalid bit depth", ERH_ERROR);

	int nb = bit_depth / 8;
	uint8_t* zdata = *p_data = calloc(nb, z_xdim * z_ydim);
	if (zdata == NULL)
		return ERROR("couldn't allocate zoomed image buffer", ERH_ERROR_ALLOC);

	for (j=0;j<z_ydim;j++)
	{
		for (i=0;i<z_xdim;i++)
		{
			int sj = j / zoom_factor;
			int si = i / zoom_factor;

			memcpy(&zdata[(j * z_xdim + i) * nb], &src[(sj * xdim + si) * nb], nb);
		}
	}

	return ERH_NO_ERROR;
}

static int get_zoom_data(int zoom_factor, int bit_depth, int* p_xdim, int* p_ydim, uint8_t** p_img_data, uint8_t** p_zoom_data)
{
	*p_zoom_data = NULL;
	if (zoom_factor == 1)			return ERH_NO_ERROR;
	if (!IN_RANGE(zoom_factor, 1, 10))	return ERROR("zoom factor must be in range [1, 10]", ERH_ERROR);

	int xdim = *p_xdim;
	int ydim = *p_ydim;
	int ret = zoom_image(bit_depth, xdim, ydim, zoom_factor, *p_img_data, p_zoom_data);
	if (ret != ERH_NO_ERROR)
		return ERROR("couldn't zoom image", ret);

	*p_xdim = xdim * zoom_factor;
	*p_ydim = ydim * zoom_factor;
	*p_img_data = *p_zoom_data;
	return ERH_NO_ERROR;
}

int write_png_binary(char* path, int xdim, int ydim, uint8_t* img_data, int zoom_factor)
{
	uint8_t* zoom_data = NULL;
	int ret = get_zoom_data(zoom_factor, 8, &xdim, &ydim, &img_data, &zoom_data);
	if (ret != ERH_NO_ERROR)
		return ERROR("couldn't get zoom data", ret);

	int bytes_per_row = NWORDS(xdim, sizeof(uint8_t));
	uint8_t* binarized = calloc(sizeof(uint8_t), bytes_per_row * ydim);
	if (binarized == NULL)
		CLEANUP("couldn't allocate binarized data", ERH_ERROR_ALLOC);

	ret = pack_bits(xdim, ydim, bytes_per_row, img_data, binarized);
	if (ret != ERH_NO_ERROR)
		CLEANUP("couldn't bitpack image", ret);

	ret = write_png_file(path, binarized, bytes_per_row, xdim, ydim, 1, PNG_COLOR_TYPE_GRAY);
	if (ret != ERH_NO_ERROR)
		CLEANUP("couldn't write PNG file", ret);

cleanup:
	free(binarized);
	free(zoom_data);
	return ret;
}

int write_png_greyscale16(char* path, int xdim, int ydim, uint16_t* gs_data, int zoom_factor)
{
	uint8_t* img_data = (uint8_t*)gs_data;
	uint8_t* zoom_data = NULL;
	int ret = get_zoom_data(zoom_factor, 16, &xdim, &ydim, &img_data, &zoom_data);
	if (ret != ERH_NO_ERROR)
		return ERROR("couldn't get zoom data", ret);

	/*int i = 0, nhits = 0;
	for (i=0;i<xdim*ydim;i++)
		nhits += img_data[i] > 0 ? 1 : 0;
	printf("nhits: %d\n", nhits);*/

	ret = write_png_file(path, img_data, sizeof(uint16_t) * xdim, xdim, ydim, 16, PNG_COLOR_TYPE_GRAY);
	if (ret != ERH_NO_ERROR)
		CLEANUP("couldn't write PNG file", ret);

cleanup:
	free(zoom_data);
	return ERH_NO_ERROR;
}

int write_png_greyscale8(char* path, int xdim, int ydim, uint8_t* img_data, int zoom_factor)
{
	uint8_t* zoom_data = NULL;
	int ret = get_zoom_data(zoom_factor, 8, &xdim, &ydim, &img_data, &zoom_data);
	if (ret != ERH_NO_ERROR)
		return ERROR("couldn't get zoom data", ret);

	ret = write_png_file(path, (uint8_t*)img_data, sizeof(uint8_t) * xdim, xdim, ydim, 8, PNG_COLOR_TYPE_GRAY);
	if (ret != ERH_NO_ERROR)
		CLEANUP("couldn't write PNG file", ret);

cleanup:
	free(zoom_data);
	return ERH_NO_ERROR;
}

int write_png_rgb(char* path, int xdim, int ydim, RGB* rgb_data, int zoom_factor)
{
	uint8_t* img_data = (uint8_t*)rgb_data;
	uint8_t* zoom_data = NULL;
	int ret = get_zoom_data(zoom_factor, 24, &xdim, &ydim, &img_data, &zoom_data);
	if (ret != ERH_NO_ERROR)
		return ERROR("couldn't get zoom data", ret);

	ret = write_png_file(path, img_data, sizeof(RGB) * xdim, xdim, ydim, 8, PNG_COLOR_TYPE_RGB);
	if (ret != ERH_NO_ERROR)
		CLEANUP("couldn't write PNG file", ret);

cleanup:
	free(zoom_data);
	return ERH_NO_ERROR;
}

