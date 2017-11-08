#ifndef UTILS_H
#define UTILS_H

//----------------------------------------------
//    includes
//----------------------------------------------
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <time.h>

//----------------------------------------------
//    defines
//----------------------------------------------
#define UNUSED(x) (void)(x)

#define KB 1024
#define MB (1024 * KB)

#define MIN(x, y) (x < y ? x : y)
#define MAX(x, y) (x > y ? x : y)

#define NDIV(n, divsize) (n / divsize + (n % divsize > 0 ? 1 : 0))
#define NWORDS(n, wordsize) NDIV(n, (wordsize * 8))

#define SIGN(x) (x >= 0 ? 1 : -1)

#define IN_RANGE(x, min, max) (x >= min && x <= max)
#define STRICTLY_IN_RANGE(x, min, max) (x > min && x < max)

#define RADIANS(x) (2.0 * M_PI * (x) / 360.0)
#define DEGREES(x) (360 * (x) / (2.0 * M_PI))

//----------------------------------------------
//    function prototypes
//----------------------------------------------
int read_file(const char* path, size_t min_filesize, uint8_t** p_buf, size_t* p_fsize);
int read_end_of_file(const char* path, size_t nbytes_to_read, uint8_t** p_buf);
int get_filesize(const char* path, size_t* p_fsize);
int dump_buffer(const char* path, size_t nbytes, uint8_t* buf);
int dump_buffer_stepped(const char* path, size_t nbytes, uint8_t* buf);
uint64_t count_set_bits(uint8_t* buf, int len);
char* copy_string(char* src);
char* path_join(char* folder, char* file);
char* str_join(char* str1, char* str2);

void setb(uint8_t* buf, int i);
int getb(uint8_t* buf, int i);
int popcount(uint8_t x);
char* filename_from_path(char* filepath);
#endif

