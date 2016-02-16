CC_MONO=g++ -Wall -mcmodel=medium

CC=$(CC_MONO)

LIBS= -lnetcdf_c++
RM=rm -rf

all: fsle2d

readparameters.o: readparameters.cpp readparameters.h
	$(CC) -c readparameters.cpp 

vectorXY.o: vectorXY.cpp vectorXY.h
	$(CC) -c vectorXY.cpp 

constants.o: constants.cpp constants.h
	$(CC) -c constants.cpp 

gridconstruction.o: gridconstruction.cpp gridconstruction.h vectorXY.o
	$(CC) -c gridconstruction.cpp

integration.o: integration.cpp integration.h vectorXY.o constants.o
	$(CC) -c integration.cpp

vflow.o: vflow.cpp vflow.h vectorXY.o constants.o
	$(CC) -c vflow.cpp -lnetcdf_c++

fsle2d.o: fsle2d.cpp
	$(CC) -c fsle2d.cpp 

fsle2d: fsle2d.o readparameters.o gridconstruction.o vflow.o constants.o vectorXY.o integration.o
	$(CC)  fsle2d.o readparameters.o gridconstruction.o vflow.o integration.o constants.o vectorXY.o -o fsle2d -lnetcdf_c++
clean:
	$(RM) *.o 
