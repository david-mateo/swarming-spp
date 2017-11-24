#include "community.h"
#include "grid.h"
#include "random.h"

/* inlines */
inline int modulo(int a, int b) {
    const int result = a % b;
    return result < 0 ? result+b: result ;
}

inline double fmodulo(double a, double b) {
    const double result = fmod(a,b);
    return result < 0. ? result+b: result ;
}


/*----------------------- Community class --------------------------*/

Community::Community(int nags , double L, Agent* ags , double* p, double* v){
    num_agents = nags ;
    agents = ags ;
    pos = p ;
    vel = v ;
    box_size = L ;
    use_grid = false ;
    grid = NULL ;
}

double* Community::get_pos(){ return pos ; }

double* Community::get_vel(){ return vel ; }

Agent* Community::get_agents(){ return agents ;}

int Community::get_num_agents(){ return num_agents ; }

double Community::get_box_size(){ return box_size ;}

// Initialization

void Community::randomize_positions(){
    for(int i=0; i<DIM*num_agents; i++)
        pos[i] = spp_random_uniform() * box_size ;
}

void Community::regular_positions(){
    int i;
    int carry ;
    int n = ceil(pow(num_agents, 1./DIM)) ; // number of agents per dimension
    int grid_indx[DIM] ;
    for(i=0; i<DIM; i++)
        grid_indx[i] = 0 ;

    for(int ia=0; ia<num_agents; ia++){
        carry = 1 ;
        for(i=0; i<DIM; i++){
            pos[ia*DIM + i] = (box_size * grid_indx[i]) / n ;
            grid_indx[i] += carry ;
            if(grid_indx[i]==n){
                grid_indx[i] = 0 ;
                carry = 1 ;
            }else{
                carry = 0 ;
            }
        }
    }
}


void Community::randomize_directions(double v0){
    /*
     * Fill the vel list with random
     * velocities of norm v0.
     */
    int i ;
    for(i=0; i<num_agents; i++)
        spp_random_vector(vel + i*DIM, v0) ;
}


// Kinematic

void Community::move(double dt){
    for(int i=0; i<num_agents*DIM; i++)
        pos[i] += dt * vel[i] ;
}

void Community::periodic_move(double dt){
    /* Range [0:box_size] in each direction */
    for(int i=0; i<num_agents*DIM; i++)
        pos[i] = fmodulo( pos[i] + dt * vel[i] , box_size );
}

// Consensus protocol

void Community::sense_velocities(double* vel_sensed){
    /*
     * If using grid, this fills the grid from scratch
     * at every iteration.
     */
    int num_neis ;
    Agent* neis ;
    if(use_grid){
        fill_grid() ;
        for(int i=0; i<num_agents; i++){
            neis = grid->get_neighborhood(agents+i , &num_neis ) ;
            agents[i].sense_velocity(num_neis , neis , vel_sensed + i*DIM) ;
        }
    }else{
        for(int i=0; i<num_agents; i++)
            agents[i].sense_velocity(num_agents , agents , vel_sensed + i*DIM) ;
    }
}

void Community::sense_noisy_velocities(double* vel_sensed){
    /*
     * If using grid, this fills the grid from scratch
     * at every iteration.
     */
    int num_neis ;
    Agent* neis  ;
    if (use_grid){
        fill_grid() ;
        for(int i=0; i<num_agents; i++){
            neis = grid->get_neighborhood(agents+i , &num_neis ) ;
            agents[i].sense_noisy_velocity(num_neis , neis , vel_sensed + i*DIM) ;
        }
    }else{
        for(int i=0; i<num_agents; i++)
            agents[i].sense_noisy_velocity(num_agents , agents , vel_sensed + i*DIM) ;
    }
}

void Community::update_velocities(double* vel_sensed){
    for(int i=0; i<num_agents*DIM; i++)
        vel[i] = vel_sensed[i] ;
}

// Printing

void Community::print_posvel(){
    int i, ia ;
    for(ia=0; ia<num_agents ; ia++){
        for(i=0; i<DIM; i++)
            printf("%f\t", pos[ia*DIM + i]) ;
        for(i=0; i<DIM; i++)
            printf("%f\t", vel[ia*DIM + i]) ;
        printf("\n") ;
    }
    printf("\n\n") ;
}

// Statistical properties

void Community::mean_position(double* meanpos){
    int i, ia ;
    for(i=0; i<DIM; i++){
        meanpos[i] = 0.0 ;
        for(ia=0; ia<num_agents; ia++){
            meanpos[i] += pos[ia*DIM + i] ;
        }
        meanpos[i] /= num_agents ;
    }
}

void Community::mean_periodic_position(double* meanpos){
    /*
     * Mean position for a periodic system
     * computed following:
     * https://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions
     */
    double mcos[DIM] , msin[DIM] ;
    int i, ia ;
    for(i=0; i<DIM; i++){
        mcos[i] = 0.0 ;
        msin[i] = 0.0 ;
    }
    for(ia=0; ia<num_agents; ia++){
        for(i=0; i<DIM; i++){
            mcos[i] += cos( pos[ia*DIM + i ] * 2. * M_PI / box_size ) ;
            msin[i] += sin( pos[ia*DIM + i ] * 2. * M_PI / box_size ) ;
        }
    }
    for(i=0; i<DIM; i++)
        meanpos[i] = box_size * ( atan2(-msin[i],-mcos[i]) + M_PI ) / ( 2. * M_PI ) ;
}

double Community::mean_velocity(double* meanvel){
    int i, ia ;
    double speed2 = 0. ;

    for(i=0; i<DIM; i++)
        meanvel[i] = 0.0 ;
    for(ia=0; ia<num_agents; ia++){
        for(i=0; i<DIM; i++)
            meanvel[i] += vel[ia*DIM + i ] ;
    }
    for(i=0; i<DIM; i++){
        meanvel[i] /= num_agents ;
        speed2 += meanvel[i] * meanvel[i] ;
    }
    return speed2 ;
}

double Community::order_parameter(double v0){
    double meanvel[DIM] ;
    double mv2 = this->mean_velocity(meanvel) ;
    return sqrt(mv2)/v0 ;
}

void Community::velocity_fluctuations(double* fluctuations, double v0){
    int i,ia ;
    double mv[DIM] ;
    double speed2, norm ;

    speed2 = this->mean_velocity(mv) ;
    norm = 1.0 / sqrt( v0 * v0 - speed2 ) ;
    for(ia=0; ia<num_agents; ia++){
        for(i=0; i<DIM; i++){
            fluctuations[i + ia * DIM] = (vel[i + ia * DIM] - mv[i]) * norm ;
        }
    }

}

void Community::correlation_histo(int n_bins, double v0, double* totalcorr, int* count){
    /*
     * Compute the correlations in velocity fluctuations in the system.
     *
     * The binsize is n_bins / (this->max_distance() * 1.000001).
     * The max_distance is slighlty increased by the factor 1.000001 to
     * avoid overflowing in
     *      bin = int( dist * bindist ) ;
     * when two agents are at _exactly_ *max_distance* from each other
     * (which can happen, for example, if they are placed in a regular grid)
     * this would be bin = n_bins but one needs bin <= n_bins-1.
     *
     * Note that the normalizing factor *norm* assumes that all the agents
     * have velocity with modulus *v0*.
     *
     */
    int i,ia,ja , bin ;
    double mv[DIM] ;
    double *v1, *v2 ;
    double dist , speed2 , norm ;
    double bindist = n_bins / (this->max_distance() * 1.000001) ;

    for(i=0; i<n_bins; i++)
        totalcorr[i] = 0.0 ;
    for(i=0; i<n_bins; i++)
        count[i] = 0 ;

    speed2 = this->mean_velocity(mv) ;
    norm = 1.0 / ( v0 * v0 - speed2 ) ;

    for(ia=0; ia<num_agents; ia++){
        v1 = agents[ia].get_vel() ;
        for(ja=ia+1; ja<num_agents; ja++){
            v2 = agents[ja].get_vel() ;
            dist = sqrt( agents[ia].distance2( agents[ja].get_pos() ) ) ;
            bin = int( dist * bindist ) ;
            count[bin] += 1 ;
                for(i=0; i<DIM; i++)
                    totalcorr[bin] += (v1[i]-mv[i]) * (v2[i]-mv[i]) ;
        }
    }
    for(i=0; i<n_bins; i++)
        totalcorr[i] *= norm ;
}


// Other
double Community::max_distance(){
    /*
     * The furthest distance in an N-dimensional
     * periodic hypercube is that between the center
     * of the cube and any of its vertices.
     * For a cube with sizes L_i, that is:
     *      d^2 = sum_i (L_i/2.)^2
     *  for L_i=L forall i,
     *      d^2 = L * N / 2^2
     */
    return box_size * sqrt(DIM / 4.) ;
}

int Community::build_network(int* num_neis, Agent*** network){
    /*
     * Note that this allocates memory to place the neighbors.
     * This is an exception to the principle that all "big" arrays
     * have to be explicitly allocated by the user, and is justified
     * because without knowing how many neighbors each agent has one
     * would need to allocate a num_agents x num_agents network.
     * Too much wasted space.
     *
     * This method DOES NOT use a Grid. The Grid deep-copies agents
     * and can screw the pointers here (especially for print_network).
     *
     */
    int jn ;
    int total_nneis = 0;
    for(int i=0; i<num_agents; i++){
        num_neis[i] = agents[i].get_neighbors(num_agents, agents) ;
        total_nneis += num_neis[i] ;
        network[i] = spp_community_alloc_neighbors(num_neis[i]) ;
        for(jn=0; jn<num_neis[i]; jn++)
            network[i][jn] = agents[i].get_neis()[jn] ;
    }
    return total_nneis ;
}

void Community::print_network(int* num_neis, Agent*** network){
    int ia, ja ;
    for(ia=0; ia<num_agents; ia++){
        for(ja=0; ja<num_neis[ia]; ja++){
            printf("%i -- %li ;\n", ia, network[ia][ja]-agents);
        }
    }
}

// Optimization related

void Community::setup_grid(Grid *g){
    grid = g ;
    use_grid = true ;
}

void Community::fill_grid(){
    grid->fill_grid( num_agents, agents ) ;
}

/*------------------- End Community class --------------------------*/

double* spp_community_alloc_space(int num_agents){
    return new double[ num_agents * DIM ] ;
}

Agent* spp_community_alloc_agents(int num_agents){
    return new Agent[num_agents] ;
}

Agent** spp_community_alloc_neighbors(int num_agents){
    return new Agent*[num_agents] ;
}

Agent* spp_community_build_agents(int num_agents, double* pos, double* vel, Agent** neis, Behavior* behavior){
    /*
     * WARNING: all the agents share the same
     * *behavior* and the same *neis* pointer.
     * Sharing the same *neis* pointer may screw with
     * paralellization!
     */
    Agent* ags = spp_community_alloc_agents(num_agents) ;
    for(int ia=0; ia<num_agents; ia++)
       ags[ia] = Agent(pos + ia*DIM, vel + ia*DIM, neis, behavior) ;
    return ags ;
}


Community spp_community_autostart(int num_agents, double speed, double box_size, Behavior* behavior){
    double* pos  = spp_community_alloc_space(     num_agents) ;
    double* vel  = spp_community_alloc_space(     num_agents) ;
    Agent** neis = spp_community_alloc_neighbors( num_agents) ;
    Agent* ags = spp_community_build_agents(num_agents, pos, vel, neis, behavior) ;
    Community com = Community(num_agents, box_size , ags, pos, vel) ;

    /* Starting positions and velocities */
    com.randomize_positions() ;
    com.randomize_directions(speed) ;
    return com ;
}
