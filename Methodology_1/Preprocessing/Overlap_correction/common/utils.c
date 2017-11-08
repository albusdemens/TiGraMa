#include "utils.h"
#include "error_handling.h"

int get_filesize(const char* path, size_t* p_fsize)
{
	*p_fsize = 0;

	FILE* fin = fopen(path, "rb");
	if (fin == NULL)
		return ERROR("couldn't open file", ERH_ERROR);

	int ret = fseek(fin, 0, SEEK_END);
	if (ret != 0)
		CLEANUP("couldn't seek to end of file", ERH_ERROR);

	*p_fsize = ftell(fin);

cleanup:
	fclose(fin);
	return ret;
}

static int read_buffer(const char* path, size_t fsize, size_t nbytes_to_read, uint8_t** p_buf)
{
	int ret = ERH_NO_ERROR;
	uint8_t* buf = NULL;
	size_t nread = 0;


	FILE* fin = fopen(path, "rb");
	if (fin == NULL)
		CLEANUP("couldn't open file", ERH_ERROR);

	ret = fseek(fin, fsize - nbytes_to_read, SEEK_SET);
	if (ret != 0)
		CLEANUP("couldn't seek to file offset", ERH_ERROR);

	buf = (uint8_t*)malloc(nbytes_to_read);
	if (buf == NULL)
		CLEANUP("couldn't allocate read buffer", ERH_ERROR_ALLOC);

	nread = fread(buf, 1, nbytes_to_read, fin);
	if (nread != nbytes_to_read)
		CLEANUP("file read failed", ERH_ERROR_ALLOC);

cleanup:
	if (ret != ERH_NO_ERROR)
	{
		free(buf);
		buf = NULL;
	}

	*p_buf = buf;
	fclose(fin);
	return ret;
}

int read_file(const char* path, size_t min_filesize, uint8_t** p_buf, size_t* p_fsize)
{
	size_t fsize = *p_fsize = 0;
	int ret = get_filesize(path, &fsize);
	if (ret != ERH_NO_ERROR)		return ERROR("couldn't get file size", ret);
	if (fsize < min_filesize)		return ERROR("file size less than minimum specified file size", ERH_ERROR);
	if (fsize > 1024 * MB)			return ERROR("file size greater than 1024 MB", ERH_ERROR);

	ret = read_buffer(path, fsize, fsize, p_buf);
	if (ret != ERH_NO_ERROR)
		return ERROR("couldn't read file", ret);

	*p_fsize = fsize;
	return ret;
}

int read_end_of_file(const char* path, size_t nbytes_to_read, uint8_t** p_buf)
{
	size_t fsize = 0;
	int ret = get_filesize(path, &fsize);
	if (ret != ERH_NO_ERROR)	return ERROR("couldn't get file size", ret);
	if (fsize < nbytes_to_read)	return ERROR("number of bytes to read is greater than filesize", ERH_ERROR);

	ret = read_buffer(path, fsize, nbytes_to_read, p_buf);
	if (ret != ERH_NO_ERROR)
		return ERROR("couldn't read file", ret);

	return ret;
}

int dump_buffer(const char* path, size_t nbytes, uint8_t* buf)
{
	FILE* fout = fopen(path, "wb");
	if (fout == NULL)
		return ERROR("couldn't open file for writing", ERH_ERROR);

	int ret = ERH_NO_ERROR;
	size_t nwritten = fwrite((char*)buf, 1, nbytes, fout);
	if (nwritten != nbytes)
		CLEANUP("whole buffer was not dumped", ERH_ERROR);

cleanup:
	fclose(fout);
	return ret;
}

int dump_buffer_stepped(const char* path, size_t nbytes, uint8_t* buf)
{
	int ret = ERH_NO_ERROR;

	FILE* fout = fopen(path, "wb");
	if (fout == NULL)
		return ERROR("couldn't open file for writing", ERH_ERROR);

	const int n = 50;
	size_t chunk_size = nbytes / n;
	size_t pos = 0;

	printf("\t__________________________________________________\n\t");

	while (nbytes)
	{
		size_t towrite = MIN(chunk_size, nbytes);

		size_t nwritten = fwrite((char*)&buf[pos], 1, towrite, fout);
		if (nwritten != towrite)
			CLEANUP("whole buffer was not dumped", ERH_ERROR);

		nbytes -= nwritten;
		pos += nwritten;
		printf("="); fflush(stdout);
	}

	printf("\n\tdone\n\n");

cleanup:
	fclose(fout);
	return ret;
}

uint64_t count_set_bits(uint8_t* buf, int len)
{
	int i = 0, j = 0;
	uint64_t count = 0;

	for (i=0;i<len;i++)
	{
		uint8_t rr = buf[i];
		for (j=0;j<8;j++)
		{
			count += rr & 1;
			rr >>= 1;
		}
	}

	return count;
}

char* copy_string(char* src)
{
	if (src == NULL)
		return NULL;

	int len = strlen(src);
	if (len == 0)
		return NULL;

	char* copy = (char*)calloc(1, len + 1);
	if (copy == NULL)
		return NULL;

	strncpy(copy, src, len);
	return copy;
}

char* path_join(char* folder, char* file)
{
	if (folder == NULL || file == NULL)
		return NULL;

	int folder_len = strlen(folder);
	int file_len = strlen(file);
	if (folder_len == 0 || file_len == 0)
		return NULL;

	char* joined = (char*)calloc(1, folder_len + file_len + 2);
	if (joined == NULL)
		return NULL;

	strncpy(joined, folder, folder_len);
	if (folder[folder_len - 1] != '/')
		joined[folder_len++] = '/';

	strncpy(&joined[folder_len], file, file_len);
	return joined;
}

char* str_join(char* str1, char* str2)
{
	if (str1 == NULL || str2 == NULL)
		return NULL;

	int len1 = strlen(str1);
	int len2 = strlen(str2);
	if (len1 == 0 || len2 == 0)
		return NULL;

	char* joined = (char*)calloc(1, len1 + len2 + 1);
	if (joined == NULL)
		return NULL;

	strncpy(joined, str1, len1);
	strncpy(&joined[len1], str2, len2);
	return joined;
}

//----------------------------------------------
//    bit-level functions
//----------------------------------------------
void setb(uint8_t* buf, int i)
{
	int wordpos = i / 8;
	int bitpos = i % 8;

	buf[wordpos] |= 1 << bitpos;
}

int getb(uint8_t* buf, int i)
{
	int wordpos = i / 8;
	int bitpos = i % 8;
	return (buf[wordpos] >> bitpos) & 1;
}

int popcount(uint8_t x)
{
	int c = 0;
	for (; x > 0; x &= x - 1)
		c++;
	return c;
}

char* filename_from_path(char* filepath)
{
	int len = strlen(filepath);
	char* filename = &filepath[len - 1];
	while (filename > filepath && *filename != '/')
		filename--;

	if (*filename == '/')
		filename++;

	return filename;
}

