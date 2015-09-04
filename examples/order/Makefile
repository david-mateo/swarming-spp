COMP= g++
#Choose library for 2D or 3D
LFLAGS= -Wall -O3 -I../../src/ -L../../src/ -lspp3d -fopenmp 
LFLAGS= -Wall -O3 -I../../src/ -L../../src/ -lspp2d -fopenmp 

vicsek_topo:	vicsek_topo.cpp
	$(COMP) -DNOISE=$(eta) $^ -o $@ $(LFLAGS) 

vicsek_metric:	vicsek_metric.cpp
	$(COMP) -DNOISE=$(eta) $^ -o $@ $(LFLAGS) 

