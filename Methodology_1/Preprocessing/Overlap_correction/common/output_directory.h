#ifndef OUTPUT_DIRECTORY_H
#define OUTPUT_DIRECTORY_H

//----------------------------------------------
//    function prototypes
//----------------------------------------------
int folder_exists(char* folder);

int mkdir_output_folder(char* output_folder);

char* create_output_folder_structure(char* output_folder, int omega_index, int alpha_index);

int get_dump_path(	char* filename_prefix, char* filename_suffix, char* output_folder,
			int alpha_index, int beta_index, int omega_index, int rindex, char** p_dumppath);

#endif

