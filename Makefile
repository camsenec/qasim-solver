.PHONY: clean

FC = gfortran
FFLAGS = -shared

.SUFFIXES: .f90 .o

.f90.o:
	$(FC)  -c $<

OBJS = constants.o \
       calc_energ.o \
			 field.o \
			 ut.o \
       qa.o

main: pimc
	time ./pimc

pimc: $(OBJS)
	$(FC) $(OBJS) -o pimc

lib: $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o libfort.so

ut.o: constants.o
calc_energ.o: constants.o ut.o
field.o: constants.o calc_energ.o ut.o
qa.o: field.o constants.o calc_energ.o ut.o

clean:
	rm -rf $(OBJS) pimc *.mod
