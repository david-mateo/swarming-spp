ROOT=/usr/local/
ROOT=$(HOME)/bin/
LIB=libspp
LIBS= $(LIB)2d.a $(LIB)3d.a
SRCS=	agent.cpp	interaction.cpp	behavior.cpp grid.cpp community.cpp \
		hostile_environment.cpp
OBJS2D=$(SRCS:.cpp=_2d.o)
OBJS3D=$(SRCS:.cpp=_3d.o)
HDRS=$(SRCS:.cpp=.h)
COMP= g++
CFLAGS= -c -Wall -O3 -ffast-math -fopenmp
LFLAGS= -Wall -O3 -ffast-math -fopenmp
CFLAGS= -c -Wall -O3 -ffast-math
LFLAGS= -Wall -O3 -ffast-math
#For debugging
#CFLAGS= -c -Wall -g  
#LFLAGS= -Wall -g

all:	$(LIBS) $(LIB).h

$(LIB)2d.a:	$(OBJS2D)
	ar -crs $@ $(OBJS2D)

$(LIB)3d.a:	$(OBJS3D)
	ar -crs $@ $(OBJS3D)

$(LIB).h:	$(HDRS)
	cat $(HDRS) | awk '!/#include/' > $@ 

%_2d.o:	%.cpp
	$(COMP) -DDIM=2 $(CFLAGS) $< -o $@  

%_3d.o:	%.cpp
	$(COMP) -DDIM=3 $(CFLAGS) $< -o $@  

install: $(LIBS) $(LIB).h
	cp $(LIB).h $(ROOT)
	cp $(LIBS) $(ROOT)
	ranlib $(ROOT)$(LIB)2d.a
	ranlib $(ROOT)$(LIB)3d.a

clean:
	rm $(OBJS2D) $(OBJS3D) $(LIBS) $(LIB).h
