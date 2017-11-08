#ifndef ERROR_HANDLING_H
#define ERROR_HANDLING_H

#define PRINT_ERRORS

//----------------------------------------------
//    macros
//----------------------------------------------
#ifdef PRINT_ERRORS
#define ERROR(msg, code) print_error(__FILE__, __PRETTY_FUNCTION__, __LINE__, msg, code)
#define CLEANUP(msg, code) {ret = code; print_error(__FILE__, __PRETTY_FUNCTION__, __LINE__, msg, code); goto cleanup;}
#else
#define ERROR(msg, code) code
#define CLEANUP(msg, code) {ret = code; goto cleanup;}
#endif

//----------------------------------------------
//    error codes
//----------------------------------------------
#define ERH_NO_ERROR	0
#define ERH_ERROR	-1
#define ERH_ERROR_ALLOC	-2
#define ERH_ERROR_PARAM	-3

//----------------------------------------------
//    function prototypes
//----------------------------------------------
int print_error(const char* file, const char* function, int line, const char* msg, int error_code);

#endif

