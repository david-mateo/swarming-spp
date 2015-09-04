#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <libspp.h>

#define NAG         5000
#define NITER       10001
#define TRANSIENT    1000
#define OUTPUT         10

#define DELTAT      1.0
#define RADIUS      1.0
#define SPEED       0.05
#define DENSITY     4.
#define BOX_SIZE    sqrt( NAG / DENSITY )
// if 3d:
// #define  BOX_SIZE    pow( NAG / DENSITY , 1./3.)

int main(int argc, char* argv[]){
    int iter ;
    double* v2    = spp_community_alloc_space( NAG) ;

    /* Set the random seed */
    long int seed ;
    if(argc>1){ seed = atol( argv[1]); }else{ seed = time(NULL) ; }
    spp_set_seed( seed ) ;

    /* Define behavior of agents */
    CartesianPeriodic g = CartesianPeriodic( BOX_SIZE ) ;
    Metric interaction = Metric( RADIUS , &g ) ;
    Vicsek_consensus behavior = Vicsek_consensus(&interaction, SPEED, NOISE) ;

    /* Create community */
    Community com = spp_community_autostart( NAG , SPEED, BOX_SIZE, &behavior ) ;

    /* Printout comments */
    printf("# Number of agents  %i\n# Metric radius     %f\n# Speed             %f\n# Noise             %f\n# Time step         %f\n# Box size          %f\n# Random seed       %li\n\n", NAG, RADIUS, SPEED, NOISE, DELTAT, BOX_SIZE, seed) ;

    /* Run some iterations to pass the
     * transient state.
     */
    for(iter=0; iter< TRANSIENT; iter++){
        com.periodic_move( DELTAT) ;
        com.sense_noisy_velocities(v2) ;
        com.update_velocities(v2) ;
    }

    /* MAIN LOOP */
    for(iter=0; iter< NITER; iter++){
        if( iter % OUTPUT == 0 ){
            printf("#Iteration: %i\tOrderpar: %f\n",iter,com.order_parameter(SPEED)) ;
            //com.print_posvel() ;
        }
        com.periodic_move( DELTAT) ;
        com.sense_noisy_velocities(v2) ;
        com.update_velocities(v2) ;
    }
    return 0;
}
