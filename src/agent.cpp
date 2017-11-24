#include "agent.h"
#include "behavior.h"
#include <stdlib.h>

Agent::Agent(double* p , double* v, Agent** ns, Behavior* bb){
    pos = p ;
    vel = v ;
    neis = ns ;
    beh = bb ;
}

// Kinetic stuff

void Agent::move(double dt){
    for(int i=0; i<DIM ; i++)
        pos[i] += vel[i]*dt ;
}

void Agent::update_vel(double* new_vel){
    /* copy values, different from vel = new_vel */
    for(int i=0; i<DIM ; i++)
        vel[i] = new_vel[i] ;
}

double* Agent::get_pos(){
    return pos ;
}

double* Agent::get_vel(){
    return vel ;
}

Agent** Agent::get_neis(){
    return neis ;
}

Behavior* Agent::get_behavior(){
    return beh ;
}

void Agent::set_pos(double* p){
    pos = p ;
}

void Agent::set_vel(double* v){
    vel = v ;
}

void Agent::set_behavior(Behavior* bb){
    beh = bb ;
}

void Agent::copy(Agent* ag){
    /* Deep copy the values of
     * ag into self.
     * TODO is *beh = *(ag->beh) always safe?
     * What if it copies from and to the
     * same memory space?
     */
    for(int i=0; i<DIM ; i++){
        pos[i] = ag->get_pos()[i] ;
        vel[i] = ag->get_vel()[i] ;
    }
    neis = ag->get_neis() ;
    *beh = *(ag->beh) ;
}

// Consensus protocol

double Agent::distance2(double *point){
    return beh->inter->g->distance2( this->pos , point) ;
}

int Agent::is_neighbor(Agent* nei){
    return beh->inter->is_neighbor( this, nei) ;
}

int Agent::get_neighbors(int n_agents, Agent* ags){
    return beh->inter->get_neighbors( this, n_agents, ags, neis);
}

void Agent::look_around(int n_agents, Agent* ags){
    beh->inter->look_around( this, n_agents, ags);
}

void Agent::sense_velocity(int num_agents, Agent* ags, double* new_vel){
    beh->sense_velocity(this, num_agents, ags, new_vel) ;
}

void Agent::sense_noisy_velocity(int num_agents, Agent* ags, double* new_vel){
    beh->sense_noisy_velocity(this, num_agents, ags, new_vel) ;
}

int Agent::sense_danger(int num_threats, Agent* threats, double* new_vel){
    return beh->sense_danger(this, num_threats, threats, new_vel) ;
}

int Agent::sense_victims(int num_agents, Agent* ags){
    return beh->sense_victims(this, num_agents, ags) ;
}

int Agent::hunt(Agent* prey, double deltat){
    return beh->hunt(this, prey, deltat) ;
}

void Agent::randomize_velocity(){
    beh->randomize_velocity(this) ;
}
