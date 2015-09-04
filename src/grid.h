#include "agent.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*
 * Class to store "Verlet lists" with the coarse position
 * of the agents in boxes.
 * The aim of this is to speed up calculations by using the
 * Grid.get_neighborhood() method so that only nearby agents
 * are used to find the neighbors of an agent.
 *
 * The computation box is divided in *nslots* regular slots
 * PER DIMENSION and a collection of agents is kept
 * in *grid* to keep track of what agents are potential
 * neighbors of an agent in a particular slot.
 * Each entry of *grid* corresponds to a slot.
 * Each entry contains a list with the agents that are
 * either in the corresponding slot or any adjacent slot.
 *
 * This tecnique makes most computations to run approximately
 * as O(N) instead of O(N**2). It is commonly used in
 * molecular dynamics calculations, see e.g.
 *  http://www.cs.ubbcluj.ro/~alibal/Teaching/Simulation/1_md_allen.pdf
 *  Sec 3.4 (page 12)
 *
 * Note that, tecnically, one can only justify using this when
 * the interaction has a finite range and the size of each
 * slot *slot_size* is greater than or equal to the interaction
 * range. For metric interaction this is the case. For topological,
 * one has to choose a number of slots conservative enough so the
 * probability of an agent having a neighbor two slots away is
 * neglibigle.
 *
 */
class Grid{
    public:
        /* Construct grid and allocate all the space it needs,
         * approximately (nslots^dimension)*max_agents*sizeof(Agent).
         * If nslots <= 3 this will print an error to stderr and
         * not allocate space, but it does NOT raise an exception.
         * Inputs:
         *      ns = nslots
         *      bs = box_size
         *      max_agents = Max number of agents that can be
         *          stored in each slot. This will affect the
         *          size allocated for grid (or rather, for
         *          each grid[i]).
         */
        Grid(int ns , double bs, int max_agents) ;
        /* Store in *ind* the n-dimensional index (i,j) or (i,j,k)
         * corresponding to a given position *pos*.
         * Inputs:
         *      pos = pointer with an n-dim position
         *      ind = pointer to store the n-dim index
         */
        void grid_index(double* pos , int *ind ) ;
        /* Return the serial index corresponding to the
         * position *pos*. This value gives the index
         * of grid[] corresponding to that position.
         */
        int  serial_index(double* pos ) ;
        /* Copy the *num_agents* contained in *ags*
         * to their corresponding slots in *grid*.
         * Each agent *ag* is _copied_ to its own
         * slot *serial_index(ag.pos)* and to all
         * the adjacent ones (8 for 2D and 26 for 3D).
         */
        void fill_grid(int num_agents, Agent* ags) ;
        /* Return a list of agents containing the
         * neighborhood of *ag* and store how many
         * agents this list contain in *num_neis*.
         * This uses the information stored in the
         * last *fill_grid* call.
         */
        Agent* get_neighborhood(Agent* ag, int* num_neis ) ;
    protected:
        /* Number of slots per dimension. */
        int nslots ;
        /* Size of the computation box */
        double box_size ;
        /* List containing a list of Agents for each slot.
         * The Grid may represent a two- or three-dimesional
         * array but grid is serialized,
         *  grid[i] = list of all agents in slot *i* and adjacent.
         */
        Agent** grid ; 
        /* List with number of agents contained in each
         * entry of *grid*
         */
        int* occupation ;
} ;


