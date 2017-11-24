#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <libspp.h>

#define NAG         1024
#define NITER      5000001
#define OUTPUT       10000
#define TRANSIENT    20000

#define DELTAT      1.0
#define SPEED       0.04
#define DENSITY     1.
#define BOX_SIZE    sqrt( NAG / DENSITY )
#define NOISE       0.05

#define NBINS       200
int main(int argc, char* argv[]){
    int iter, bin ;
    double* v2    = spp_community_alloc_space(NAG ) ;
    double totalcorr[NBINS] ;
    int count[NBINS] ;
    double maxdis ;
    int num_neis[NAG] ;
    Agent** network[NAG] ;
    double mean_neis ;

    /* Set the random seed */
    long int seed ;
    if(argc>1){ seed = atol( argv[1]); }else{ seed = time(NULL) ; }
    spp_set_seed( seed ) ;


    CartesianPeriodic g = CartesianPeriodic( BOX_SIZE ) ;
    Metric interaction = Metric( RADIUS , &g) ;
    Vicsek_consensus behavior = Vicsek_consensus(&interaction, SPEED, NOISE) ;

    /* Create community */
    Community com = spp_community_autostart( NAG , SPEED, BOX_SIZE, &behavior) ;
    maxdis = com.max_distance() ;

    /* Set the agents in a grid and freeze the interaction. */
    com.regular_positions() ;
    mean_neis = com.build_network(num_neis, network) * 1.0 / NAG ;
    NetworkInteraction grid_interaction = NetworkInteraction(com.get_agents(), num_neis, network, &g) ;
    behavior.inter = &grid_interaction ;

    /* Printout comments */
    printf("# Number of agents  %i\n# Grid connections  %f\n# Speed             %f\n# Noise             %f\n# Time step         %f\n# Box size          %f\n# Random seed       %li\n\n", NAG, mean_neis, SPEED, NOISE, DELTAT, BOX_SIZE, seed) ;

    for(iter=0; iter< TRANSIENT; iter++){
        com.sense_noisy_velocities(v2) ;
        com.update_velocities(v2) ;
    }

    for(iter=0; iter< NITER; iter++){
        if( iter % OUTPUT == 0 ){
            printf("#Iteration: %i\tOrderpar: %f\n",iter,com.order_parameter(SPEED)) ;
            com.correlation_histo(NBINS, SPEED, totalcorr, count) ;
            for(bin=0; bin< NBINS; bin++)
                printf("%f\t%f\t%i\n", (bin+0.5)*maxdis/NBINS ,totalcorr[bin], count[bin]) ;
            printf("\n\n") ;
        }
        com.sense_noisy_velocities(v2) ;
        com.update_velocities(v2) ;
    }
    return 0;
}
