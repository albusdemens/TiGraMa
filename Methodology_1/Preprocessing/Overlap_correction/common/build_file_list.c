#include <fnmatch.h>
#include <fts.h>
#include <errno.h>
#include <stdio.h>
#include <stdbool.h>
#include "error_handling.h"
#include "file_list.h"
#include "output_directory.h"


static bool valid_file(FTS* tree, FTSENT* node, char* filename_filter)
{
	if (node->fts_level > 0 && node->fts_name[0] == '.')
	{
		fts_set(tree, node, FTS_SKIP);
		return false;
	}

	if (node->fts_level > 1)	//don't search recursively
		return false;

	if (!(node->fts_info & FTS_F))
		return false;

	if (filename_filter == NULL)
		return true;

	return fnmatch(filename_filter, node->fts_name, 0) == 0 ? true : false;
}

int build_file_list(FILE_LIST* flist, char* input_folder, char* filename_filter)
{
	int ret = ERH_NO_ERROR;
	int fileno = 0;
	FTSENT* node = NULL;
	char* paths[2] = {input_folder, NULL};


	if (!folder_exists(input_folder))
		return ERROR("folder doesn't exist", ERH_ERROR);

	file_list_init(flist);

	FTS* tree = fts_open(paths, FTS_NOCHDIR, 0);
	if (tree == NULL)
		CLEANUP("couldn't open folder for traversal", ERH_ERROR);

	while ((node = fts_read(tree)))
	{
		if (valid_file(tree, node, filename_filter))
		{
			ret = file_list_add_item(flist, node->fts_path, node->fts_name);
			if (ret != ERH_NO_ERROR)
				goto cleanup;

			fileno++;
			if (fileno % 10000 == 0)
				printf("files found: %d\n", fileno);
		}
	}

	printf("total files found: %d\n", flist->count);
	if (flist->count == 0)
		CLEANUP("no files matching filename filter were found", ERH_ERROR);

cleanup:
	if (fts_close(tree))
		return ERROR("fts_close error", ERH_ERROR);

	return ret;
}

