# CC 	   = g++
CC 	   = mpic++ -std=c++11 -g

#CFLAGS     = -pg -fopenmp
CFLAGS     = -O3 -Wall -I${HOME}/local/fftw_mpi -I${HOME}/local/fftw_mpi/include -L/home/chris/local/fftw_mpi -L/home/chris/local/fftw_mpi/include -I${FFTWHOME}/include

# CFLAGS     = -O3 -openmp -math -Wall -I/home/hchao/Install/fftw3_openmp/include
# LIBS      = -openmp  -lm  -O3 -lfftw3_omp -lfftw3 -lpthread -L/home/hchao/Install/fftw3_openmp/lib
LIBS      = -lm  -O3  -lfftw3 -lfftw3_mpi -lpthread -L/home/chris/local/fftw_mpi/lib -L${HOME}/local/fftw_mpi -L${HOME}/local/fftw_mpi/include -L/home/chris/local/fftw_mpi -L/home/chris/local/fftw_mpi/include  -I${FFTWHOME}/include

C         = mpic++ -std=c++11
# CC       = ~/usr/local/bin/mpic++ -std=c++11
#FFTW_LOC = /opt/seas/pkg/gcc/fftw3/mpi/double/3.3.7
FFTWHOME = ${HOME}/local/fftw_mpi
#FFTW_LOC = include2/fftw_mpi
#EIGEN_LOC = include2/eigen
EIGEN_PRINT_LOC = ${HOME}/Install/gdb
INCLUDE = -I${FFTWHOME}/include
CFLAGS     = -I${FFTW_LOC}/include -I${FFTWHOME}/include -I${EIGEN_LOC} -O3 -Wno-unused-result -Wno-write-strings
LIBS      = -lm -lfftw3_mpi -lfftw3 -O3 -lpthread  -L${FFTWHOME}/lib
DLIBS      = -lm -lfftw3_mpi -lfftw3  -L${FFTWHOME}/lib
DFLAGS     = -I${FFTW_LOC}/include  -I${FFTWHOME}/include -I${EIGEN_PRINT_LOC} -I${EIGEN_LOC} -Wno-unused-result -Wno-write-strings

#############################################################################
# nothing should be changed below here

SRCS = main.cpp array_utils.cpp die.cpp  random.cpp grid_utils.cpp \
			 fftw_wrappers.cpp initialize.cpp config_utils.cpp io_utils.cpp \
			 update_positions.cpp forces.cpp integ_utils.cpp read_input.cpp \
			 bonded.cpp calc_unb.cpp stress.cpp 
       
       
			 


OBJS = ${SRCS:.cpp=.o}

.cpp.o:
	${CC} ${CFLAGS} ${DFLAGS} -c  $<

dmft:  ${OBJS}
	$(CC) ${CFLAGS} ${DFLAGS} -o $@ ${OBJS} $(LIBS)

clean:
	rm -f *.o
	rm -f dmft
	rm -f *~

