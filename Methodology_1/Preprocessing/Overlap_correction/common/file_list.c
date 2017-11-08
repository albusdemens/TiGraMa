#include "file_list.h"
#include "error_handling.h"
#include "utils.h"


#define FILE_LIST_MAGIC 0x34256199

void file_list_init(FILE_LIST* list)
{
	list->magic = FILE_LIST_MAGIC;
	list->count = 0;
	list->capacity = 0;
	list->file = NULL;
}

int file_list_uninit(FILE_LIST* list)
{
	if (list == NULL)			return ERROR("file list argument is NULL", ERH_ERROR);
	if (list->magic != FILE_LIST_MAGIC)	return ERROR("argument is not an initialized file list", ERH_ERROR);

	int i = 0;
	for (i=0;i<list->count;i++)
	{
		free(list->file[i].path);
		free(list->file[i].name);
	}

	list->capacity = 0;
	list->count = 0;
	free(list->file);
	list->file = NULL;
	return ERH_NO_ERROR;
}

int file_list_add_item(FILE_LIST* list, char* path, char* name)
{
	if (list == NULL)			return ERROR("file list argument is NULL", ERH_ERROR);
	if (list->magic != FILE_LIST_MAGIC)	return ERROR("argument is not an initialized file list", ERH_ERROR);
	if (path == NULL)			return ERROR("path argument is NULL", ERH_ERROR);
	if (name == NULL)			return ERROR("name argument is NULL", ERH_ERROR);

	if (list->count == list->capacity)
	{
		int ncap = list->capacity == 0 ? 256 : 2 * list->capacity;

		void* ptr = realloc(list->file, sizeof(FILE_INFO) * ncap);
		if (ptr == NULL)
			return ERROR("list realloc failed", ERH_ERROR_ALLOC);

		list->file = (FILE_INFO*)ptr;
		list->capacity = ncap;
	}

	FILE_INFO* file = &list->file[list->count];

	file->path = copy_string(path);
	file->name = copy_string(name);
	if (file->path == NULL)			return ERROR("couldn't allocate filepath string buffer", ERH_ERROR_ALLOC);
	if (file->name == NULL)			return ERROR("couldn't allocate filename string buffer", ERH_ERROR_ALLOC);

	list->count++;
	return ERH_NO_ERROR;
}

