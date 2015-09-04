#include "grid.h"

#if DIM==2
#define NSLOTSD  ( nslots * nslots )
#elif DIM==3
#define NSLOTSD  ( nslots * nslots * nslots )
#endif

Grid::Grid( int ns , double bs , int max_agents){
    /* 
     * Note that if ns<=3 this prints an error and
     * doesn't allocate memory but it does not
     * raise an exception.
     * Erratic behavior can occur if called with ns<=3.
     */
    int is ;
    if( ns <= 3){
        fprintf(stderr,"WARNING Grid: Invalid number of slots %i\n", ns) ;
        return ;
    }
    nslots = ns ;
    box_size = bs ;
    occupation = new int[ NSLOTSD ] ;
    grid       = new Agent*[ NSLOTSD ] ;
    for(is=0; is< NSLOTSD ; is++){
        occupation[is] = 0 ;
        grid[is] = new Agent[ max_agents ] ;
    }
}

void Grid::grid_index(double* pos, int *ind){
    for(int i=0 ; i<DIM ; i++)
        ind[i] = floor( pos[i] / box_size * nslots ) ;
}

int Grid::serial_index(double* pos){
#if DIM==2
    return floor( pos[0] / box_size * nslots ) * nslots +
           floor( pos[1] / box_size * nslots ) ;
#elif DIM==3
    return floor( pos[0] / box_size * nslots ) * nslots * nslots +
           floor( pos[1] / box_size * nslots ) * nslots +
           floor( pos[2] / box_size * nslots ) ;
#endif
}

void Grid::fill_grid(int num_agents, Agent* agents){
    /* 
     * Put every agent in its cell and all
     * the adjacent ones, meaning all that have an
     * n-dimensioanl index with compoenents equal or
     * +-1 different from the cell.
     *
     */
    int  serial_ind, ind[DIM] ;
    for(int is=0 ; is < NSLOTSD ; is++)
        occupation[is] = 0 ;

    for(int ia=0 ; ia < num_agents ; ia++){
        grid_index( agents[ia].get_pos(), ind ) ;
#if DIM==2
        for(int i=ind[0]-1 ; i<=ind[0]+1 ; i++){
            for(int j=ind[1]-1 ; j<=ind[1]+1 ; j++){
                serial_ind = (i<0?i+nslots:i%nslots) * nslots + 
                             (j<0?j+nslots:j%nslots) ;
                grid[serial_ind][occupation[serial_ind]] = agents[ia] ; 
                occupation[serial_ind] += 1 ;
            }
        }
#elif DIM==3
        for(int i=ind[0]-1 ; i<=ind[0]+1 ; i++){
            for(int j=ind[1]-1 ; j<=ind[1]+1 ; j++){
                for(int k=ind[2]-1 ; k<=ind[2]+1 ; k++){
                    serial_ind = (i<0?i+nslots:i%nslots) * nslots * nslots +
                                 (j<0?j+nslots:j%nslots) * nslots +
                                 (k<0?k+nslots:k%nslots) ;
                    grid[serial_ind][occupation[serial_ind]] = agents[ia] ; 
                    occupation[serial_ind] += 1 ;
                }
            }
        }
#endif
    }
}

Agent* Grid::get_neighborhood(Agent* ag, int* num_neis ){
    int index = this->serial_index( ag->get_pos() ) ;
    *num_neis = occupation[index] ;
    return grid[index] ;
}
