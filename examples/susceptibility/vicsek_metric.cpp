#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <libspp.h>

#define NAG         1024
#define NITER       750000
#define OUTPUT        1000
#define TRANSIENT    20000

#define DELTAT      1.0
#define SPEED       0.04
#define DENSITY     1.
#define BOX_SIZE    sqrt( NAG / DENSITY )
#define NOISE       0.05

#define NBINS       200

int main(int argc, char* argv[]){
    int iter, bin ;
    double* v2    = spp_community_alloc_space( NAG) ;
    double totalcorr[NBINS] ;
    int count[NBINS] ;
    double maxdis ;

    /* Set the random seed */
    long int seed ;
    if(argc>1){ seed = atol( argv[1]); }else{ seed = time(NULL) ; }
    spp_set_seed( seed ) ;

    /* Printout comments */
    printf("# Number of agents  %i\n# Metric radius     %f\n# Speed             %f\n# Noise             %f\n# Time step         %f\n# Box size          %f\n# Random seed       %li\n\n", NAG, RADIUS, SPEED, NOISE, DELTAT, BOX_SIZE, seed) ;

    CartesianPeriodic g = CartesianPeriodic( BOX_SIZE ) ;
    Metric interaction = Metric( RADIUS , &g) ;
    Vicsek_consensus behavior = Vicsek_consensus(&interaction, SPEED, NOISE) ;

    /* Create community */
    Community com = spp_community_autostart( NAG , SPEED, BOX_SIZE, &behavior) ;
    maxdis = com.max_distance() ;
    /* Use grid */
    int nslots = (int) BOX_SIZE / RADIUS ;
    if (nslots > 50 )
        nslots = 50 ; // Above 50 slots there's little speed gain in the calculation.
    Grid* grid ;
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
            printf("#Iteration: %i\tOrderpar: %f\n",iter,com.order_parameter(SPEED)) ;
            com.correlation_histo(NBINS, SPEED, totalcorr, count) ;
            for(bin=0; bin< NBINS; bin++)
                printf("%f\t%f\t%i\n", (bin+0.5)*maxdis/NBINS ,totalcorr[bin], count[bin]) ;
            printf("\n\n") ;
        }
        com.periodic_move( DELTAT) ;
        com.sense_noisy_velocities(v2) ;
        com.update_velocities(v2) ;
    }
    return 0;
}
