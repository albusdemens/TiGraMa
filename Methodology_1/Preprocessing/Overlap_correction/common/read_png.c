#include "utils.h"
#include "error_handling.h"
#include <png.h>

static int read_png_file_internal(const char* input_path, int* p_xdim, int* p_ydim, uint8_t** p_img_data, int* p_bit_depth)
{
	int width, height;
	png_byte color_type;
	png_byte bit_depth;

	png_structp png_ptr;
	png_infop info_ptr;
	//int number_of_passes;

	unsigned char header[8];    // 8 is the maximum size that can be checked

	//open file and test for it being a png
	FILE* fp = fopen(input_path, "rb");
	if (fp == NULL)
		return ERROR("[read_png_file] input file could not be opened for reading", ERH_ERROR);

	size_t nread = fread(header, 1, 8, fp);
	if (nread != 8 || png_sig_cmp(header, 0, 8))
	{
		fclose(fp);
		return ERROR("[read_png_file] input file is not recognized as a PNG file", ERH_ERROR);
	}

	//initialize stuff
	png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	if (!png_ptr)
	{
		fclose(fp);
		return ERROR("[read_png_file] png_create_read_struct failed", ERH_ERROR);
	}

	info_ptr = png_create_info_struct(png_ptr);
	if (!info_ptr)
	{
		fclose(fp);
		return ERROR("[read_png_file] png_create_info_struct failed", ERH_ERROR);
	}

	//if (setjmp(png_jmpbuf(png_ptr)))
	//{
	//	fclose(fp);
	//	return ERROR("[read_png_file] Error during init_io", ERH_ERROR);
	//}

	png_init_io(png_ptr, fp);
	png_set_sig_bytes(png_ptr, 8);

	png_read_info(png_ptr, info_ptr);

	width = png_get_image_width(png_ptr, info_ptr);
	height = png_get_image_height(png_ptr, info_ptr);
	color_type = png_get_color_type(png_ptr, info_ptr);
	*p_bit_depth = bit_depth = png_get_bit_depth(png_ptr, info_ptr);

	if (color_type != 0)
	{
		fclose(fp);
		return ERROR("input file colour type != 0", ERH_ERROR);
	}

	*p_xdim = width;
	*p_ydim = height;

	//number_of_passes = png_set_interlace_handling(png_ptr);
	png_read_update_info(png_ptr, info_ptr);

	//read file
	//if (setjmp(png_jmpbuf(png_ptr)))
	//{
	//	fclose(fp);
	//	return ERROR("[read_png_file] Error during read_image", ERH_ERROR);
	//}

	//big-endian/little-endian correction
	png_set_swap(png_ptr);

	//read data
	png_bytep* row_pointers = (png_bytep*)malloc(sizeof(png_bytep) * height);
	if (row_pointers == NULL)
	{
		fclose(fp);
		return ERROR("couldn't allocate row pointers", ERH_ERROR_ALLOC);
	}

	int bytes_per_row = NWORDS(width * bit_depth, sizeof(uint8_t));
	uint8_t* img_data = *p_img_data = calloc(bytes_per_row, height);
	if (img_data == NULL)
	{
		fclose(fp);
		free(row_pointers);
		return ERROR("couldn't allocate image data buffer", ERH_ERROR_ALLOC);
	}

	int i = 0;
	for (i=0;i<height;i++)
		row_pointers[i] = (png_byte*)&img_data[i * bytes_per_row];

	png_read_image(png_ptr, row_pointers);

//cleanup:
	free(row_pointers);

	png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
	fclose(fp);
	return ERH_NO_ERROR;
}

static int getb_png(uint8_t* buf, int i)
{
	int wordpos = i / 8;
	int bitpos = 7 - i % 8;
	return (buf[wordpos] >> bitpos) & 1;
}

static void expand_bits(int xdim, int ydim, uint8_t* binary_data, uint8_t* expanded)
{
	int i = 0, j = 0;
	int bytes_per_row = NWORDS(xdim, sizeof(uint8_t));

	for (j=0;j<ydim;j++)
		for (i=0;i<xdim;i++)
			expanded[j * xdim + i] = getb_png(&binary_data[j * bytes_per_row], i);
}

int read_png_file(const char* input_path, int* p_xdim, int* p_ydim, uint8_t** p_img_data, int* p_bit_depth)
{
	int xdim = 0, ydim = 0, bit_depth = 0;
	uint8_t* img_data = NULL;

	int ret = read_png_file_internal(input_path, &xdim, &ydim, &img_data, &bit_depth);
	if (ret != ERH_NO_ERROR)
		return ERROR("couldn't read png file", ret);

	if (bit_depth == 1)
	{
		uint8_t* expanded = calloc(sizeof(uint8_t), xdim * ydim);
		if (expanded == NULL)
			CLEANUP("couldn't allocate expanded buffer", ERH_ERROR_ALLOC);

		expand_bits(xdim, ydim, img_data, expanded);
		free(img_data);
		img_data = expanded;
	}

cleanup:
	*p_xdim = xdim;
	*p_ydim = ydim;
	*p_img_data = img_data;
	*p_bit_depth = bit_depth;
	return ret;
}

