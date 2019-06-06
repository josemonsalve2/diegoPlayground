# Makefile DEMAC GeMM
PROGRAM = main
GCC = mpicc
NPP = 4

${PROGRAM}: ${PROGRAM}.cpp
	${GCC} -o ${PROGRAM}.o ${PROGRAM}.cpp
run:
	mpiexec -np ${NPP} ./${PROGRAM}.o
clean:
	rm -f ${PROGRAM}.o M_gen_* M_out_*