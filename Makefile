.PHONY: clean

FC = gfortran
FFLAGS = -shared

.SUFFIXES: .f90 .o

.f90.o:
	$(FC)  -c $<

OBJS = constants.o \
       calc_energ.o \
			 field.o \
			 scalar_pointer_char_wrapper.o\
       qa.o

main: pimc
	time ./pimc

pimc: $(OBJS)
	$(FC) $(OBJS) -o pimc

lib: $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o libfort.so

calc_energ.o: constants.o
field.o: constants.o calc_energ.o
qa.o: field.o constants.o calc_energ.o scalar_pointer_char_wrapper.o

clean:
	rm -rf $(OBJS) pimc *.mod
