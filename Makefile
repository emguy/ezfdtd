EXECUTABLE=ezfdtd
SRC = ezfdtd.c domain.c ade.c step.c excitation.c probes.c dft.c pml.c cpml.c classical.c tools.c mem.c h5io.c mur.c
HEADERS = $(SRC:.c=.h)
OBJECTS = $(SRC:.c=.o)
OBJECTS_EXE = $(SRC:.c=.wo)
CC_EXE=/usr/MinGW-w64/bin/x86_64-w64-mingw32-gcc -I/usr/MinGW-w64/include -L/usr/MinGW-w64/lib
CFLAGS=-march=core-avx-i -pipe -O2
CFLAGS_EXE=-O2

all: cscope options $(EXECUTABLE)

LIBMATH=/usr/MinGW-w64/mingw/lib/libm.a

cscope:
	@echo update cscope cache
	@cscope -Rbq

options:
	@echo fdtd_2d build options:
	@echo "CC           = $(CC)"
	@echo "CFLAGS       = ${CFLAGS}"
	@echo "CFLAGS_EXE   = ${CFLAGS_EXE}"
	@echo "LDFLAGS      = ${LDFLAGS}"
	@echo "CC_EXE       = $(CC_EXE)"

%.o: %.c %.h
	@echo $(CC) $<
	$(CC) -Wall -c  -D ADD_COLOR ${CFLAGS} $< 

%.wo: %.c %.h
	@echo $(CC_EXE) -Wall -c $<
	$(CC_EXE) -Wall -c ${CFLAGS_EXE} $< -o $@

$(EXECUTABLE): $(OBJECTS)
	@echo CC -o $@
	@$(CC) -Wall -o $@ $(OBJECTS) -lhdf5 -lm
	@cp -vf ezfdtd ../../working_dir/FDTD_test/fdtdio/ezfdtd

exe: $(OBJECTS_EXE)
	@echo CC_EXE -o $@
	@$(CC_EXE) -Wall -o $(EXECUTABLE).exe $(OBJECTS_EXE) $(LIBMATH)  -lhdf5  -lz

clean:
	@rm -f $(OBJECTS) $(OBJECTS_EXE) $(EXECUTABLE) $(EXECUTABLE).exe cscope*

.PHONY: all options clean
