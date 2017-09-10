#include "community.h"

/*
 * Extension of Community to handle
 * systems containing two kinds of agents:
 * predators and preys.
 * The *ags* (preys) and the *preds* (predator)
 * should have the appropiate behavior with
 * sense_danger() for preys, and sense_victims()
 * and hunt() for predators.
 *
 * This class is designed with a small number
 * of predators in mind. It should work with
 * any number of predators, but it does not
 * implement the same performance improvements
 * for large number of agents for the predators.
 * TODO Test this with number of predators >1.
 *
 */
class HostileEnvironment : public Community{
    public:
        /* Construct the HostileEnvironment class.
         * Inputs:
         *      nags = number of agents (preys)
         *      in the community.
         *      L = size of box. It is fixed
         *      to be the same in all directions
         *      for facilitating the use of a Grid.
         *      ags = pointer to the array storing
         *      the *nags* Agent instances (preys).
         *      p = pointer to the array to store
         *      the positions of agents
         *      (size depens on DIM).
         *      v = pointer to the array to store
         *      the velocities of agents
         *      (size depens on DIM).
         *      npreds = number of agent predators
         *      in the community.
         *      ags = pointer to the array storing
         *      the *npred* Agent instances
         *      (predators).
         */
        HostileEnvironment(int nags , double L, Agent* ags , double* p, double* v, int npreds , Agent* preds) ;
        /* Return the pointer to the array of agents (predators).*/
        Agent* get_predators() ;
        /* Same as Community::sense_velocity but also
         * calls the sense_danger() method of *ags*.
         * Returns the number of agents fleeing from
         * a predator.
         */
        int sense_velocities_danger(double* vel_sensed) ;
        /* same as sense_velocities_danger() but using
         * sense_noisy_velocity() instead of sense_velocity().
         */
        int sense_noisy_velocities_danger(double* vel_sensed) ;
        /* Move the predator towards the closest prey. If
         * the predator can catch the prey, remove the prey
         * from the community using the remove_dead() function.
         * For some application, you may want to change this
         * function so that it calls replace_dead() instead.
         * Returns the number of agents removed.
         */
        int hunt(double dt) ;
        /* Remove *ags[ia]* from the community. This function
         * reduces *num_agents* by one and rearranges the agents
         * in *ags*. The size of *ags* does not change.
         */
        void remove_dead(int ia) ;
        /* Replace *ags[ia]* by a new agent with consensus
         * velocity, generated as far as possible from
         * the dead one.
         */
        virtual void replace_dead(int ia) ;
        /* Same as Community::print_posvel() but for
         * the position and velocity of the predators.
         */
        void print_predators_posvel() ;
    protected:
        /* Number of predators. */
        int num_predators ;
        /* Array of predator agents
         *      Size: num_predators
         */
        Agent* predators ;
} ;

/* Similar as spp_community_autostart.
 * Return a HostileEnvironment instance "ready to use" from scratch.
 * Allocate the required space, initialize the agents
 * and randomize their positions and velocities.
 * The speed *speed* is only used to randomize velocities
 * with a constant norm.
 * The predators start at a random position with 0 velocity.
 *
 * Uses spp_community_alloc_neighbors and the returned Community
 * is therefore not thread-safe for parallelization.
 */
HostileEnvironment spp_hostile_autostart(int num_agents, double speed, double box_size, Behavior* agsbeh, int num_predators, Behavior* predsbeh) ;
