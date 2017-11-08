#include <sys/stat.h>
#include <errno.h>
#include "utils.h"
#include "error_handling.h"

int folder_exists(char* folder)
{
	struct stat st;
	return stat(folder, &st) == 0;
}

int mkdir_output_folder(char* output_folder)
{
	struct stat st;

	if (stat(output_folder, &st) == 0)	//folder present
		return ERH_NO_ERROR;

	if (errno != ENOENT)			//ENOENT: no such file or directory
		return ERROR("couldn't stat directory", ERH_ERROR);

	if (mkdir(output_folder, 0777) != 0 && errno != EEXIST)
		return ERROR("mkdir failed", ERH_ERROR);

	return ERH_NO_ERROR;
}

char* create_output_folder_structure(char* output_folder, int omega_index, int alpha_index)
{
	char* temp = NULL;
	char str[16] = {0};

	sprintf(str, "w=%d", omega_index);
	temp = path_join(output_folder, str);
	if (temp == NULL)
		return NULL;

	mkdir_output_folder(temp);
	free(temp);

	sprintf(str, "w=%d/a=%d", omega_index, alpha_index);
	temp = path_join(output_folder, str);
	if (temp == NULL)
		return NULL;

	mkdir_output_folder(temp);
	return temp;
}

int get_dump_path(char* filename_prefix, char* filename_suffix, char* output_folder,
			int alpha_index, int beta_index, int omega_index, int rindex, char** p_dumppath)
{
	int ret = ERH_NO_ERROR;
	char* output_dir = NULL;

	char filename[64] = {0};
	sprintf(filename, "%s_w=%d_a=%d_b=%d_i=%d%s", filename_prefix, omega_index, alpha_index, beta_index, rindex, filename_suffix);

	output_dir = create_output_folder_structure(output_folder, omega_index, alpha_index);
	if (output_dir == NULL)
		CLEANUP("couldn't create output folder", ERH_ERROR);

	*p_dumppath = path_join(output_dir, filename);
	if (*p_dumppath == NULL)
		CLEANUP("couldn't allocate dump path", ERH_ERROR_ALLOC);

cleanup:
	free(output_dir);
	return ret;
}

