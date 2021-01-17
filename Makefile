TARGET = fem
OBJECTS = linarg.o fem2d.o
MOD_FILE = linarg.mod

FC = gfortran

.SUFFIXES: .o .f90

%.o: %.f90
	${FC} -O2 -c $<

%.mod: %.f90 %.o
	@:

${TARGET}: ${OBJECTS}
	${FC} -O2 -o $@ ${OBJECTS}
