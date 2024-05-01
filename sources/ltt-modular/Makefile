SOURCE_FOLDER = $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))

EXE     = ltt-modular.exe
FC      = gfortran
FCFLAGS = -fcheck=all -fbacktrace -fmax-errors=10 -I${ECCODES_INSTALL_ROOT}/include
LDFLAGS = -L${ECCODES_INSTALL_ROOT}/lib -Wl,-rpath=${ECCODES_INSTALL_ROOT}/lib -leccodes_f90
OBJS    = ltt_main.o module_computeDelays.o module_data_types.o module_io.o module_nwpProjection.o module_skyview.o module_undulations.o module_utility.o

vpath %.f90 $(SOURCE_FOLDER)

$(EXE) : $(OBJS)
	$(FC) -o $@ $(FCFLAGS) $(LDFLAGS) $^

%.o : %.f90
	$(FC) -c $(FCFLAGS) $<

clean :
	rm -f $(EXE) *.o *.mod

ltt_main.o             : module_data_types.o module_io.o module_skyview.o module_nwpProjection.o module_computeDelays.o module_undulations.o
module_computeDelays.o : module_data_types.o module_utility.o module_io.o module_nwpProjection.o module_undulations.o module_utility.o
module_io.o            : module_data_types.o
module_nwpProjection.o : module_data_types.o module_utility.o module_undulations.o
module_skyview.o       : module_data_types.o
module_undulations.o   : module_data_types.o


