COMP= g++
LFLAGS= -Wall -O3 -I../../src/ -L../../src/ -lspp3d
LFLAGS= -Wall -O3 -I../../src/ -L../../src/ -lspp2d

vicsek_topo_k%:	vicsek_topo.cpp
	$(COMP) -DOUTDEGREE=$* $^ -o $@ $(LFLAGS)
	
vicsek_metric_r%:	vicsek_metric.cpp
	$(COMP) -DRADIUS=$* $^ -o $@ $(LFLAGS)
