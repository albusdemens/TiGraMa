# c-compiler
CC=clang++
CFLAGS=-stdlib=libstdc++ -g -O3 -Wall -Wextra
LDLIBS= -lcfitsio -lm -lpng -fopenmp

#folder where intermediate object files will be stored
OBJDIR=objdir

#plain c-code files
CDEPS=FitsFile.o TPX_Cube.o TPX_CubeRead.o utils.o output_directory.o build_file_list.o file_list.o error_handling.o

COBJ = $(addprefix $(OBJDIR)/, $(CDEPS))

all: $(OBJDIR) recon

$(OBJDIR):
	mkdir $(OBJDIR)

recon: $(OBJDIR) $(COBJ)
	$(CC) $(CFLAGS) $(COBJ) /usr/local/Cellar/cfitsio/3.370/lib/libcfitsio.a -o fits_correction $(LDLIBS)

############################
#    main
############################
$(OBJDIR)/FitsFile.o: FitsFile.cpp
	$(CC) -c $(CFLAGS) FitsFile.cpp $(LDLIBS) -o $@

$(OBJDIR)/TPX_Cube.o: TPX_Cube.cpp
	$(CC) -c $(CFLAGS) TPX_Cube.cpp $(LDLIBS) -o $@

$(OBJDIR)/TPX_CubeRead.o: TPX_CubeRead.cpp
	$(CC) -c $(CFLAGS) TPX_CubeRead.cpp $(LDLIBS) -o $@

############################
#    common
############################
$(OBJDIR)/utils.o: ../common/utils.c
	$(CC) -c $(CFLAGS) ../common/utils.c $(LDLIBS) -o $@

$(OBJDIR)/output_directory.o: ../common/output_directory.c
	$(CC) -c $(CFLAGS) ../common/output_directory.c $(LDLIBS) -o $@

$(OBJDIR)/build_file_list.o: ../common/build_file_list.c
	$(CC) -c $(CFLAGS) ../common/build_file_list.c $(LDLIBS) -o $@

$(OBJDIR)/file_list.o: ../common/file_list.c
	$(CC) -c $(CFLAGS) ../common/file_list.c $(LDLIBS) -o $@

$(OBJDIR)/math_utils.o: ../common/math_utils.c
	$(CC) -c $(CFLAGS) ../common/math_utils.c $(LDLIBS) -o $@

$(OBJDIR)/error_handling.o: ../common/error_handling.c $(OBJDIR)/utils.o
	$(CC) -c $(CFLAGS) ../common/error_handling.c $(LDLIBS) -o $@

############################
#    clean
############################
clean:
	rm -rf $(OBJDIR) fits_correction

