#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <libspp.h>

#define NAG         2048
#define NITER       5000
#define TRANSIENT   2000
#define OUTPUT       100

#define DELTAT      0.2
#define NOISE       0.05
#define DRAD        0.5
#define SPEED       0.2 
#define VPRED       1.5 * SPEED 
#define DENSITY     1.
#define BOX_SIZE    sqrt( NAG / DENSITY )


int main(int argc, char* argv[]){
    int iter ;
    double* v2    = spp_community_alloc_space( NAG) ;
    double* dist2 = spp_community_alloc_space( NAG) ;
    int death, avoidance_time = 0;

    /* Set the random seed */
    long int seed ;
    if(argc>1){ seed = atol( argv[1]); }else{ fprintf(stderr,"No seed specified.\n\tUsage: ./vicsek_metric seed\n"); abort() ;}
    spp_set_seed( seed ) ;

    /* Printout comments */
    printf("# Number of agents  %i\n# Outdegree         %i\n# Speed             %f\n# Noise             %f\n# Time step         %f\n# Box size          %f\n# Random seed       %li\n\n", NAG, OUTDEGREE, SPEED, NOISE, DELTAT, BOX_SIZE, seed) ;

    /* Define behavior of agents */
    CartesianPeriodic g = CartesianPeriodic( BOX_SIZE ) ;
    Topologic interaction = Topologic( OUTDEGREE , &g , dist2) ;
    // Vicsek_prey* behavior = new Vicsek_prey[NAG] ;
    Vicsek_prey     prey_beh = Vicsek_prey( &interaction, SPEED, NOISE, DRAD );
    Vicsek_predator pred_beh = Vicsek_predator( &interaction, VPRED, NOISE );

    /* Create community */
    HostileEnvironment com = spp_hostile_autostart( NAG , SPEED, BOX_SIZE, &prey_beh , 1 , &pred_beh ) ;
    Grid* grid ;
    /* Use grid */
    int nslots = (int) sqrt( (0.5 * NAG ) / OUTDEGREE ) ; 
    if( nslots > 3 ){
        try{
            grid = new Grid( nslots , BOX_SIZE , NAG ) ;
            com.setup_grid( grid) ;
            printf("# Using grid with %i slots/dim.\n", nslots) ;
        }catch(...){
        printf("# Failed to use grid with %i slots/dim.\n", nslots) ;
        }
    }else{
        printf("# Too few slots per dimension (%i).\n", nslots) ;
    }
 
    for(iter=0; iter< TRANSIENT; iter++){
        com.periodic_move( DELTAT) ;
        com.sense_noisy_velocities(v2) ;
        com.update_velocities(v2) ;
    }
 
    for(iter=0; iter< NITER; iter++){
        if( iter % OUTPUT == 0 ){
            printf("#Iteration: %i\tNum agents: %i\n", iter, com.get_num_agents() ) ;
        }
        death = com.hunt( DELTAT) ;
        if(death){
            printf("%i\n", avoidance_time) ;
            avoidance_time = 0 ;
        }else{
            avoidance_time +=1 ;
        }
        com.periodic_move( DELTAT) ;
        com.sense_noisy_velocities_danger(v2) ;
        com.update_velocities(v2) ;
    }
    return 0;
}
