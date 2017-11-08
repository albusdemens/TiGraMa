#ifndef FILE_LIST_H
#define FILE_LIST_H

typedef struct
{
	char* path;
	char* name;
} FILE_INFO;

typedef struct
{
	int magic;
	int count;
	int capacity;
	FILE_INFO* file;
} FILE_LIST;

//----------------------------------------------
//    function prototypes
//----------------------------------------------
void file_list_init(FILE_LIST* list);
int file_list_uninit(FILE_LIST* list);
int file_list_add_item(FILE_LIST* list, char* path, char* name);

#endif

