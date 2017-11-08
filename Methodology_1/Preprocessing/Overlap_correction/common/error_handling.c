#include <stdio.h>
#include <stdlib.h>

int print_error(const char* file, const char* function, int line, const char* msg, int error_code)
{
	printf("\n\nerror\tfile: %s\n", file);
	printf("\tline: %d\n", line);
	printf("\tfunction: %s\n", function);
	printf("\terror message: %s\n", msg);
	printf("\terror code: %d\n", error_code);
	(void)fflush(stdout);

//exit(error_code);
	return error_code;
}

