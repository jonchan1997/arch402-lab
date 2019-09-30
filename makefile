## Compiler, tools and options
CC         = gcc
LINKER     = gfortran
CCFLAGS    = -O1 -mavx -O3 #-Q --help=optimizers
LDFLAGS    = 

LAPACKHOME = /home/youngjon/lapack-3.8.0
INCPATH    =  -I$(LAPACKHOME)/CBLAS/include
LIBS       =  -L$(LAPACKHOME) -lcblas -lrefblas 


## Files
OBJECTS = c_timer.o lab.o 
TARGET  = lab


## Implicit rules
.SUFFIXES: .c
.c.o:
	$(CC) -c  $(CCFLAGS) $(INCPATH) $<


## Build rules
all: $(TARGET)


$(TARGET): $(OBJECTS)
	$(LINKER) -o $@  $(OBJECTS)  $(CCFLAGS) $(LDFLAGS) $(LIBS)

clean:
	rm -f $(OBJECTS) $(TARGET)
	rm -f *~ core
