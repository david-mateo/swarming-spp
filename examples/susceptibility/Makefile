COMP= g++
LFLAGS= -Wall -O3 -I../../src/ -L../../src/ -lspp3d -ffast-math
LFLAGS= -Wall -O3 -I../../src/ -L../../src/ -lspp2d -ffast-math

vicsek_topo_k%:	vicsek_topo.cpp
	$(COMP) -DOUTDEGREE=$* $^ -o $@ $(LFLAGS)
	
vicsek_metric_r%:	vicsek_metric.cpp
	$(COMP) -DRADIUS=$* $^ -o $@ $(LFLAGS)

vicsek_regular_r%:	vicsek_regular.cpp
	$(COMP) -DRADIUS=$* $^ -o $@ $(LFLAGS)
