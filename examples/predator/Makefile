COMP= g++
LFLAGS= -Wall -O3 -I../../src/ -L../../src/ -lspp3d -ffast-math
LFLAGS= -Wall -O3 -I../../src/ -L../../src/ -lspp2d -ffast-math

predator_topo_k%:	predator_topo.cpp
	$(COMP) -DOUTDEGREE=$* $^ -o $@ $(LFLAGS)
	
predator_metric_r%:	predator_metric.cpp
	$(COMP) -DRADIUS=$* $^ -o $@ $(LFLAGS)
