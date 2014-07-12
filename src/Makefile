# Comment with icc
CC=gcc
CFLAGS= -O3 -Wall -Wshadow -pedantic -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -DVER32 -fopenmp \
	-I/opt/local/include/ -L/opt/local/lib/   -g  #-p  #-m64

# Comment with gcc
#CC=icc
#CFLAGS= -O3 -fast -Wall -Wshadow -pedantic -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -DVER32 \
	-I/opt/local/include/ -L/opt/local/lib/   #-g  #-p  #-m64

# The source files, object files, libraries and executable name.
SRCFILES= mlRho.c eprintf.c stringUtil.c interface.c mlComp.c piComp.c profile.c ld.c profileTree.c deltaComp.c filterPro.c
OBJFILES= mlRho.o eprintf.o stringUtil.o interface.o mlComp.o piComp.o profile.o ld.o profileTree.o deltaComp.o filterPro.o
LIBS= -lgsl -lgslcblas -lm

EXECFILE= mlRho
DIRECTORY= MlRho

VERSION= 2.7

# The make rule for the executable
.PHONY : all
all : $(EXECFILE)

$(EXECFILE) : $(OBJFILES)
	$(CC) $(CFLAGS) -fopenmp -o $(EXECFILE) $(OBJFILES) $(LIBS)

# Other Standard make rules
lint : 
	lint $(SRCFILES) | more

clean:
	rm -f *.o *~ mlRho

