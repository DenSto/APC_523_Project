CC = icc
MPICC = mpic++
MPIFLAGS = -DUSE_MPI=0

CFLAGS= -g  -Wall 
CFLAGS= -O3 -ipo -xhost -qopt-report=5 -qopt-report-phase=vec -xCORE-AVX2
#CFLAGS= -O3

include Makefile.obj

SRCDIR = $(PWD)

LCONFIG_DIR = $(SRCDIR)/../libconfig-1.5
LCONFIG_INCL = $(LCONFIG_DIR)/include
LCONFIG_LIB = $(LCONFIG_DIR)/lib/libconfig++.a

LIBS = $(LCONFIG_LIB)
LIBS_PAR = $(LCONFIG_LIB)
LIBS_PIC = $(LCONFIG_LIB)

PROGRAM = EMOOPIC
PROGRAM_PIC = EMOOPIC_PIC
PROGRAM_PAR = $(PROGRAM)_par

all: $(PROGRAM) $(PROGRAM_PIC)

$(PROGRAM): $(OBJ) $(LIBS)
	$(CC) -o $@ $^ $(CFLAGS) -lm -lz -ldl 

$(PROGRAM_PIC): $(OBJ_PIC) $(LIBS_PIC)
	$(CC) -o $@ $^ $(CFLAGS) -lm -lz -ldl -DPART_IN_CELL

$(PROGRAM_PAR): $(OBJ_PAR) $(LIBS_PAR)
	$(MPICC) -o $@ $^ $(CFLAGS) -lm -lz -ldl $(MPIFLAGS)

%.o: %.cpp $(LIBS)
	$(CC) -c -o $@ $< $(CFLAGS)

%_par.o : %.cpp $(LIBS_PAR)
	$(MPICC) -c -o $@ $< $(CFLAGS) $(MPIFLAGS) 

%_pic.o : %.cpp $(LIBS_PIC)
	$(CC) -c -o $@ $< $(CFLAGS) -DPART_IN_CELL

$(LCONFIG_LIB):
	cd ..; \
	tar xf libconfig-1.5.tar.gz; \
	cd libconfig-1.5; \
	./configure --prefix=$(LCONFIG_DIR); \
	make; \
	make install; \
	cd $(SRCDIR)

clean:
	rm -f $(PROGRAM) $(PROGRAM_PAR) $(PROGRAM_PIC) $(OBJ) $(OBJ_PAR) $(OBJ_PIC) .depend

fullclean: clean
	rm -rf ../libconfig-1.5
run_tests:
	cd ../test; \
	make; \
	make run_tests; \
	cd ../src;

