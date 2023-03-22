PROGS   = main
OBJECTS = $(PROGS).o fftw_f.o NSE.o file_management.o
HEADERS = $(PROGS).h fftw_f.h NSE.h file_management.h



IDIR =include
CC=mpicc
CFLAGS=-I -Wall -g -I/usr/local/include -L/usr/local/lib

ODIR =obj


LIBS= -lfftw3l -lm 

_DEPS = $(HEADERS)
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = $(OBJECTS)
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(PROGS): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f *.o *~ core $(INCDIR)/*~           \
         rm -f $(IDIR)/*.o *~ core $(INCDIR)/*~   \
         rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~    \
         rm -f $(PROGS)


IDIT =analysis/files

deletefiles:
	rm -f *.o *~ core $(INCDIR)/*~           \
		 rm -f $(IDIT)/*.txt *~ core $(INCDIR)/*~   \
		 rm -f $(PROGS)

PLT =analysis/plots
cleanall:
	rm -f *.o *~ core $(INCDIR)/*~           \
         rm -f $(IDIR)/*.o *~ core $(INCDIR)/*~   \
         rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~    \
         rm -f $(PROGS)\
         rm -f $(IDIT)/*.txt *~ core $(INCDIR)/*~   \
         rm -f $(PLT)/*.gif *~ core $(INCDIR)/*~   \

