#include "hostile_environment.h"
#include "grid.h"

/*----------------------- Hostile class --------------------------*/
HostileEnvironment::HostileEnvironment(int nags , double L, Agent* ags , double* p, double* v, int npreds , Agent* preds) : Community(nags, L, ags, p, v) {
    num_predators = npreds ;
    predators = preds ;
}

Agent* HostileEnvironment::get_predators(){
    return predators ;
}

int HostileEnvironment::sense_velocities_danger(double* vel_sensed){
    int ia , fleeing = 0;

    this->sense_velocities(vel_sensed) ;

    for(ia=0; ia<num_agents; ia++)
        fleeing += agents[ia].sense_danger(num_predators, predators, vel_sensed + ia*DIM) ;
    return fleeing ;
}

int HostileEnvironment::sense_noisy_velocities_danger(double* vel_sensed){
    int ia , fleeing = 0 ;

    this->sense_noisy_velocities(vel_sensed) ;

    for(ia=0; ia<num_agents; ia++)
        fleeing += agents[ia].sense_danger(num_predators, predators, vel_sensed + ia*DIM) ;
    return fleeing ;
}

int HostileEnvironment::hunt(double dt){
    int ip , iprey , deaths = 0 , kill;
    for(ip=0; ip<num_predators; ip++){
        iprey = predators[ip].sense_victims(num_agents, agents) ;
        kill  = predators[ip].hunt( agents+iprey , dt) ;
        if(kill){
            deaths +=1 ;
            this->remove_dead(iprey) ;
            // this->replace_dead(iprey) ;
        }
    }
    return deaths ;
}

void HostileEnvironment::remove_dead(int ia){
    /* To avoid messing with dynamic memory allocation,
     * this function simply reduces the number of
     * agents *num_agents* and copies the last agent
     * stored in *ags[num_agents]* into the "dead"
     * agent *ags[ia]*. Note the copy is a deep copy
     * of the state of the agent, not of the reference.
     * This avoid messing up the correpondence of *pos*
     * and *vel* to *ags*.
     */
    num_agents -= 1 ;
    agents[ia].copy( agents + num_agents )  ;
}

void HostileEnvironment::replace_dead(int ia){
    /* Replace the agent ia by a new agent placed at
     * the opposite end of the periodic space and
     * with random velocity.
     */
    int i ;
    /* displace the agent to the opposite end of the box */
    for(i=ia*DIM; i<(ia+1)*DIM; i++)
        pos[i] = fmodulo( pos[i] + 0.5*box_size , box_size );
    /* set the velocity to the sensed velocity in the new location */
    agents[ia].randomize_velocity() ;
}

void HostileEnvironment::print_predators_posvel(){
    int i, ip ;
    for(ip=0; ip<num_predators ; ip++){
        for(i=0; i<DIM; i++)
            printf("%f\t", predators[ip].get_pos()[i]) ;
        for(i=0; i<DIM; i++)
            printf("%f\t", predators[ip].get_vel()[i]) ;
        printf("\n") ;
    }
}


/*------------------- End Hostile class --------------------------*/

HostileEnvironment spp_hostile_autostart(int num_agents, double speed, double box_size, Behavior* agsbeh, int num_predators, Behavior* predsbeh){
    /* Agents (preys) */
    double* pos  = spp_community_alloc_space(     num_agents) ;
    double* vel  = spp_community_alloc_space(     num_agents) ;
    Agent** neis = spp_community_alloc_neighbors( num_agents) ;
    Agent* ags = spp_community_build_agents(num_agents, pos, vel, neis, agsbeh) ;
    /* Predators */
    double* ppos = spp_community_alloc_space( num_predators ) ;
    double* pvel = spp_community_alloc_space( num_predators ) ;
    Agent* preds = spp_community_build_agents(num_predators, ppos, pvel, neis, predsbeh) ;

    HostileEnvironment hos = HostileEnvironment(num_agents, box_size , ags, pos, vel, num_predators, preds) ;

    for(int i=0; i< num_predators * DIM; i++){
       ppos[i] = spp_frandom() * box_size ;
       pvel[i] = 0. ;
    }
    /* Starting positions and velocities */
    hos.randomize_positions() ;
    hos.randomize_directions(speed) ;
    return hos ;
}
