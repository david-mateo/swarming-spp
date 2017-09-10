#include "behavior.h"
#include <stdlib.h>
#include <stdio.h>

class Grid ;

/*
 * Community class implemented to easily
 * manage a collection of Agent instances.
 * On one hand, this provides a cleaner
 * interface for the user to interact with
 * the ensamble of agents by calling the
 * methods of Community instead of looping
 * for each agent and computing global
 * properties such as correlations or order
 * parameter.
 * On the other hand, this provides some
 * performance improvements by making sure
 * all the positions and velocities are
 * aligned in a contiguos memory segment and
 * using molecular dynamics-like tricks to
 * speed up calculations by the use of a Grid
 * so that only nearby agents are checked for
 * neighborhood.
 *
 * The Community class manages several arrays
 * of data but it does not allocate any space
 * of order *nags* (some small arrays are
 * allocated for temporary results).
 * The user is responsible for all the space
 * allocation, and can use the functions provided
 * here to automatize this process.
 *
 */
class Community{
    public:
        /* Construct the community class.
         * Inputs:
         *      nags = number of agents in
         *      the community.
         *      L = size of box. It is fixed
         *      to be the same in all directions
         *      for facilitating the use of a Grid.
         *      ags = pointer to the array storing
         *      the *nags* Agent instances.
         *      p = pointer to the array to store
         *      the positions of agents
         *      (size depens on DIM).
         *      v = pointer to the array to store
         *      the velocities of agents
         *      (size depens on DIM).
         */
        Community(int nags , double L, Agent* ags , double* p, double* v) ;
        /* Return the pointer to the position array.
         * The position of agent *i* corresponds to
         * the values
         *      get_pos()[i*DIM : i*DIM + DIM]
         */
        double* get_pos() ;
        /* Return the pointer to the velocity array.
         * The velocity of agent *i* corresponds to
         * the values
         *      get_pos()[i*DIM : i*DIM + DIM]
         */
        double* get_vel() ;
        /* Return the pointer to the array of agents.*/
        Agent* get_agents() ;
        /* Return how many agents are in the community.*/
        int get_num_agents() ;
        /* Return box size.*/
        double get_box_size() ;
        /* Store random values in the position array *pos*.
         * The values for each dimension are in the
         * [0:box_size] range.
         */
        void randomize_positions();
        /* Store random values in the velocity array *vel*.
         * Each random velocity has fixed norm *v0*.
         */
        void randomize_directions(double v0) ;
        /* Move all the agents during *dt* time, meaning
         * pos += dt * vel
         */
        void move(double dt) ;
        /* Move all the agents during *dt* time, assuming
         * periodic boundary conditions in each direction.
         */
        void periodic_move(double dt) ;
        /* Sense the consensus velocities using
         * each agent's behavior and store the
         * result in *vel_sensed*.
         * Calls Agent->behavior->sense_velocity
         * and uses a Grid if setup. If using
         * Grid, each call to this also fills the grid.
         * Note: it is not safe to use the class'
         * own *vel* as *vel_sensed*.
         */
        void sense_velocities(double* vel_sensed) ;
        /* Same as sense_velocities but calls the
         * Agent->behavior->sense_nosiy_velocity
         * method.
         * Note: it is not safe to use the class'
         * own *vel* as *vel_sensed*.
         */
        void sense_noisy_velocities(double* vel_sensed) ;
        /* Copy the values in *vel_sensed* to *vel*.
         * This needs to be done separate from the sense_*
         * method to make sure the velocities are
         * updated synchronously.
         */
        void update_velocities(double* vel_sensed) ;
        /* Print the position and velocity of each
         * agent to stdin. The format used is
         * x    y   vz  vz                  (2D)
         * x    y   z   vz  vy  vz          (3D)
         * Prints two blank lines at the end.
         */
        void print_posvel() ;
        /* Store the mean position (center of mass)
         * of all the agents in *meanpos*.
         */
        void mean_position(double* meanpos) ;
        /* Store the mean position (center of mass)
         * of all the agents in *meanpos*. Takes
         * into account the periodic boundary
         * conditions to compute a "mean angular position".
         * See https://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions
         */
        void mean_periodic_position(double* meanpos) ;
        /* Store the mean velocity
         * of all the agents in *meanvel*.
         * Return the norm^2 of the mean.
         */
        double mean_velocity(double* meanvel) ;
        /* Return the order parameter of velocity alignment,
         * defined as
         *      op = sqrt{|(1/N)\sum_{i=0}^N v_i|^2} / v0
         *  (the norm of the mean velocity divided by *v0*.)
         *  If the norm of each velocity is =< v0, then the
         *  order parameter is a number between 0 and 1.
         *  The higher the order parameter, the more aligned
         *  the agents are.
         */
        double order_parameter(double v0) ;
        /* Compute the correlations in velocity fluctuations in the system.
         * Mathematical formulation based on the work by Attanasi et al. in
         *      PLoS Comput Biol 10, e1003697 (2014)
         *
         * The function computes the total correlation *totalcorr* and the
         * count of pairs of agents at a certain distance *count* separate
         * instead of returning just *totalcorr/count* (Eq 2) to be able
         * to compute correctly the cumulative correlation (Eq 3) and
         * the susceptibility.
         * The total correlation at the bin *i*, *totalcorr[i]*, corresponds
         * to the sum of correlations between agents at a distance *r* such
         * that
         *      i * max_dist() / n_bins < r < (i+1) * max_dist() / n_bins
         * The *count[i]* is the number of agents at such a distance from
         * each other.
         *
         * WARNING: The normalizing factor *norm* assumes that all the agents
         * have velocity with modulus *v0*.
         */
        void correlation_histo(int n_bins, double v0, double* totalcorr, int* count) ;
        /* Return the distance between the two farthest points in
        * the computation box with periodic boundary conditions.
        */
        double max_distance() ;
        /* Start using a Grid to compute
         * the neighbors of each agent.
         * The Grid instance *g* has to
         * be initialized by the user before
         * calling this function.
         */
        void setup_grid(Grid* g) ;
        /* Fill the grid calling its own
         * fill_grid method. This method
         * does NOT check that if grid has
         * been setup or not.
         */
        void fill_grid() ;
    protected:
        /* Number of agents. */
        int num_agents ;
        /* positions of the agents
         *      Size: num_agents * DIM
         */
        double* pos ;
        /* velocities of the agents
         *      Size: num_agents * DIM
         */
        double* vel ;
        /* Array of agents
         *      Size: num_agents
         */
        Agent* agents ;
        /* Box size, same in all directions. */
        double box_size ;
        /* False by default. Turns to True
         * when a Grid instance is inseted
         * in Community through setup_grid().
         */
        bool use_grid ;
        /* Grid instance to store the coarse
         * location of each agent. See Grid
         * documentation for more info.
         * Is it initialized to NULL and can
         * be set with setup_grid().
         */
        Grid* grid ;
} ;

// Utils for automatization of the setup of a Community.

/* Allocate space for storing a vector per agent,
 * i.e. num_agents * DIM doubles, and return the
 * pointer to the array.
 */
double* spp_community_alloc_space(int num_agents) ;
/* Allocate space for a vector of Agent instances
 * and return the pointer to the array.
 */
Agent* spp_community_alloc_agents(int num_agents) ;
/* Allocate space for a vector of Agent POINTERS
 * and return the pointer to the array.
 */
Agent** spp_community_alloc_neighbors(int num_agents) ;
/* Construct an array of *num_agents* Agents where
 * each agent uses consequetive slots of pos (vel) to
 * store its position (velocity). All the agents
 * have the same behavior *behavior* and the same
 * list to store pointers to neighbors *neis*.
 *
 * WARNING: All the agents sharing the same *neis*
 * can yield to problems if paralelization is used.
 *
 * TODO add a spp_community_build_agents_parallel(...)
 */
Agent* spp_community_build_agents(int num_agents, double* pos, double* vel, Agent** neis, Behavior* behavior) ;
/* Return a Community instance "ready to use" from scratch.
 * Allocate the required space, initialize the agents
 * and randomize their positions and velocities.
 * The speed *speed* is only used to randomize velocities
 * with a constant norm.
 *
 * Uses spp_community_alloc_neighbors and the returned Community
 * is therefore not thread-safe for parallelization.
 */
Community spp_community_autostart(int num_agents, double speed, double box_size, Behavior* behavior) ;


/* inlines */
inline int modulo(int a, int b) {
    const int result = a % b;
    return result < 0 ? result+b: result ;
}

inline double fmodulo(double a, double b) {
    const double result = fmod(a,b);
    return result < 0. ? result+b: result ;
}
