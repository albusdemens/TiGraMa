#ifndef READ_PNG_H
#define READ_PNG_H

//----------------------------------------------
//    function prototypes
//----------------------------------------------
int read_png_file(const char* input_path, int* p_xdim, int* p_ydim, uint8_t** p_img_data, int* p_bit_depth);

#endif

