class Agent ;
#include "interaction.h"

/*
 * Abstract Behavior class used as a template
 * for different Behavior implementations.
 * Every implementation defines a particular
 * consensus protocol, i.e. a rule to determine
 * what velocity the agent should align to given
 * a set of connected agents or 'neighbors'.
 * Every behavior must implement at least:
 *
 *      sense_velocity: determine the new velocity
 *          *new_vel* that agent *ag* should align
 *          to given *num_agents* potential neighbors
 *          stored in *ags*. Which of them are
 *          actual neighbors is computed via *inter*.
 *      sense_noisy_velocity: same as sense_velocity
 *          with an added random noise.
 *
 * The other functions are optional depending on
 * the desired behavior.
 *
 */
class Behavior {
    public:
        /* pure virtual, must be implemented */
        virtual void sense_velocity(Agent* ag, int num_agents, Agent* ags, double* new_vel) = 0;
        /* pure virtual, must be implemented */
        virtual void sense_noisy_velocity(Agent* ag, int num_agents, Agent* ags, double* new_vel) = 0;
        /* pure virtual, must be implemented */
        virtual void randomize_velocity(Agent* ag) = 0;
        /* optional */
        virtual int sense_danger(Agent* ag, int num_threats, Agent* threats, double* new_vel) {return 0;};
        /* optional */
        virtual int sense_victims(Agent* ag, int num_agents, Agent* ags) {return 0;} ;
        /* optional */
        virtual int hunt(Agent* ag, Agent* prey, double deltat) {return 0;};
        /* Interaction pointer that determines if a
         * given agent is a neighbor of another given
         * agent. See documentation of Interaction
         * for more info.
         */
        Interaction* inter ;
    protected:
        /* Noise level parameter used by sense_noisy_velocity.
         * Assumed to be in the [0:1] range.
         */
        double noise ;
} ;

/*
 * Original Vicsek consensus protocol
 * from Vicsek et al. PRL 75, 1226 (1995).
 * The consensus velocity is the mean velocity
 * of all the neighbors scaled to have a
 * constant norm v0. Noise is introduced as
 * a random rotation on the computed new
 * velocity.
 * The agents in the original model explicitly
 * include themselves in the mean. Here, the
 * decision of whether an agent is a neighbor
 * of itself or not is left to the Interaction
 * instance *inter*.
 * (it is, for Metric and Topologic)
 *
 */
class Vicsek_consensus : public Behavior {
    public:
        Vicsek_consensus() {} ;
        Vicsek_consensus(Interaction* ii, double v0, double noise) ;
        /* Store the mean velocity of *ag*'s neighbors in *new_vel*,
         * re-scaled to have a *v0* norm.
         */
        void sense_velocity(Agent* ag, int num_agents, Agent* ags, double* new_vel) ;
        /* Rotate a vector *v* by a random angle between [-noise*pi : noise*pi].
         * For dimensions higher than 2 a random rotation axis is also chosen, which
         * increases considerably the amount of computation required for this.
         */
        void rotate(double* v) ;
        /* Sense velocity using sense_velocity() and then rotate the
         * sensed velocity *new_vel* using the rotate() method.
         */
        void sense_noisy_velocity(Agent* ag, int num_agents, Agent* ags, double* new_vel) ;
        /* Sets the velocity of *ag* to a random vector with norm v0. */
        void randomize_velocity(Agent* ag) ;
    protected:
        /* Fixed norm of the agent velocity.
         */
        double v0 ;
} ;

/*
 * Same as Vicsek consensus but the
 * noise is implemented as a vector added
 * to the consensus velocity instead of an
 * angle rotating the velocity.
 * The model is discussed in
 *      Chate et al. Phys. Rev. E 77, 046113  (2008)
 */
class Chate_consensus : public Vicsek_consensus {
    public:
        Chate_consensus() {} ;
        Chate_consensus(Interaction* ii, double v0, double noise) : Vicsek_consensus(ii, v0, noise) {} ;
        /* Sense velocity including a vectorial noise term
         * proportional to the number of neighbors following
         * the receipe from:
         *      PHYSICAL REVIEW E 77, 046113  2008
         */
        void sense_noisy_velocity(Agent* ag, int num_agents, Agent* ags, double* new_vel) ;
} ;

/*
 * Vicsek consensus protocol plus
 * a prey behavior where the agents
 * sense danger and if detected act
 * on it, ignoring the consensus
 * protocol.
 * A threat is detected if any of
 * the *num_threats* agents in *threats*
 * is at a distance2 less than
 * *dectection_radius2*. If detected,
 * the agent proceds to align itself
 * in order to radially flee from
 * the first threat detected.
 *
 */
class Vicsek_prey : public Vicsek_consensus {
    public:
        Vicsek_prey() : Vicsek_consensus() {} ;
        Vicsek_prey(Interaction* ii, double v0, double noise, double dradius) ;
        /* If a threat is detected, return 1 and store the flee velocity
         * in *new_vel*. Else, do nothing and return 0.
         */
        int sense_danger(Agent* ag, int num_threats, Agent* threats, double* new_vel) ;
    protected:
        /* Square of the maximum distance at which
         * agents are capable of detecting threats.
         *
         * TODO For the sake of generality,
         * this should be replaced by an
         * Interaction instance that determines
         * if the agent detects the threat or not.
         */
        double detection_radius2 ;
} ;

/*
 * Vicsek consensus protocol plus a predator
 * behavior where the agent senses victims
 * and hunts an agent if it close enough.
 *
 *      hunt: if *prey* is closer than *v0*x*deltat*
 *          to the predator *ag*, then *ag*
 *          moves to *prey* position and the
 *          function returns 1. Otherwise,
 *          move *ag* for *v0*x*deltat* in
 *          *prey* direction and return 0.
 *          Note the actual "hunting" is left
 *          for the user or the HostileEnvironment
 *          class to do.
 *
 */
class Vicsek_predator : public Vicsek_consensus {
    public:
        Vicsek_predator() : Vicsek_consensus() {} ;
        Vicsek_predator(Interaction* ii, double v0, double noise) : Vicsek_consensus(ii, v0, noise) {} ;
        /* Find which agent out of the *num_agents*
         * stored in *ags* is closer to the
         * predator *pred*.
         * Return the index of that agent.
         */
        int sense_victims(Agent* pred, int num_agents, Agent* ags) ;
        /* If *pred* can move to *prey*'s position
         * in a time *deltat* or less, then *pred*
         * moves to where *prey* is and return 1.
         * Else, move *pred* in the direction
         * of *prey* and return 0.
         * Return 1 if the succesfully hunted
         * and 0 otherwise.
         */
        int hunt(Agent* pred, Agent* prey, double deltat) ;
} ;
