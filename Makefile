CPP=g++
CPPFLAGS=-O2 -Wall

LDFLAGS= -lnetcdf_c++ 
RM=rm -rf

DESTDIR=$(HOME)/bin

common_src=  constants.cpp vflow.cpp vectorXY.cpp vectorIJ.cpp rparameters.cpp gridconstruction.cpp constants.cpp integration.cpp
common_obj=$(common_src:.cpp=.o) 
common_dep=$(common_obj:.o=.d)  # one dependency file for each source

#RTIME2D 
rtime2d_src= rtime2d.cpp
rtime2d_obj=$(rtime2d_src:.cpp=.o) 
rtime2d_dep=$(rtime2d_obj:.o=.d)  # one dependency file for each source

#FSLE2D
fsle2d_src= fsle2d.cpp
fsle2d_obj=$(fsle2d_src:.cpp=.o) 
fsle2d_dep=$(fsle2d_obj:.o=.d)  # one dependency file for each source

.PHONY: all rtime2d fsle2d

all: rtime2d fsle2d

rtime2d: $(common_obj) $(rtime2d_obj)
	$(CPP) -o $@ $^ $(LDFLAGS)

fsle2d: $(common_obj) $(fsle2d_obj)
	$(CPP) -o $@ $^ $(LDFLAGS)

-include $(common_dep) # include all dep files in makefile
-include $(rtime2d_dep) 
-include $(fsel2d_dep) 


%.d: %.cpp 	# rule to generate a dep file by using the g++ prepocesor
	$(CPP) $(CPPFLAGS) -MM -MT $(@:.d=.o) $< -MF $@

%.o: %.cpp
	$(CPP) $(CPPFLAGS) -o $@ -c $<

.PHONY: debug rtime2d fsle2d
debug: CPPFLAGS+= -DDEBUG -ggdb # debug with gdb
debug: rtime2d fsle2d

.PHONY: clean
clean:
	$(RM) $(common_obj) $(rtime2d_obj) $(fsle2d_obj) *.d *~ *# rtime2d fsle2d

.PHONY: install
install: rtime2d fsle2d
	install $^ $(DESTDIR)
