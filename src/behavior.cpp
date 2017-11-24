#include "behavior.h"
#include "agent.h"
#include "random.h"

/*
 * Vicsek Consensus
 */
Vicsek_consensus::Vicsek_consensus(Interaction* ii, double vzero, double ns) {
    inter = ii ;
    v0 = vzero ;
    noise = ns ;
}

void Vicsek_consensus::sense_velocity(Agent* ag, int num_agents, Agent* ags, double* new_vel){
    int i, j ;
    double  v2 = 0.0 ;
    int num_neis ;
    Agent** neis = ag->get_neis() ;

    for(i=0; i<DIM ; i++)
        new_vel[i] = 0.;
    num_neis = inter->get_neighbors(ag, num_agents, ags, neis) ;

    for(j=0; j<num_neis ; j++){
        for(i=0; i<DIM ; i++) new_vel[i] += neis[j]->get_vel()[i];
    }

    for(i=0; i<DIM ; i++) v2 += new_vel[i]*new_vel[i] ;
    for(i=0; i<DIM ; i++) new_vel[i] *= v0/sqrt(v2) ;
}

void Vicsek_consensus::rotate( double *v){
    double theta = noise * 2.0 * M_PI * (spp_random_uniform()-0.5) ;
#if DIM==2
    double tmp ;
    tmp  = cos(theta) * v[0] - sin(theta) * v[1] ;
    v[1] = sin(theta) * v[0] + cos(theta) * v[1] ;
    v[0] = tmp ;
#elif DIM>2
    double axis[DIM] ;
    double av = 0.0; /* = axis*v */
    int i ;
    spp_random_vector(axis, 1.0) ; // unitary vector
    for(i=0; i<DIM ; i++) av += v[i]*axis[i] ;

    for(i=0; i<DIM ; i++){
        v[i] = cos(theta)*v[i] + sin(theta)/sqrt( v0*v0 - av*av)*( axis[i]*v0*v0 - v[i]*av) ;
    }
#endif
}

void Vicsek_consensus::sense_noisy_velocity(Agent* ag, int num_agents, Agent* ags, double* new_vel){
    this->sense_velocity(ag, num_agents, ags, new_vel) ;
    this->rotate(new_vel) ;
}

void Vicsek_consensus::randomize_velocity(Agent* ag){
    spp_random_vector(ag->get_vel(), v0) ;
}


/*
 * Chate Consensus
 * (same as Vicsek but with vectorial noise)
 */
void Chate_consensus::sense_noisy_velocity(Agent* ag, int num_agents, Agent* ags, double* new_vel){
    int i, j ;
    double  v2 = 0.0 ;
    int num_neis ;
    Agent** neis = ag->get_neis() ;
    num_neis = inter->get_neighbors(ag, num_agents, ags, neis) ;

    /* Start the vel to a random vector
     * of norm noise*v0*num_neis */
    spp_random_vector(new_vel, noise*v0*num_neis) ;

    for(j=0; j<num_neis ; j++){
        for(i=0; i<DIM ; i++) new_vel[i] += neis[j]->get_vel()[i];
    }

    for(i=0; i<DIM ; i++) v2 += new_vel[i]*new_vel[i] ;
    for(i=0; i<DIM ; i++) new_vel[i] *= v0/sqrt(v2) ;
}


/*
 * Vicsek Prey
 */
Vicsek_prey::Vicsek_prey(Interaction* ii, double vzero, double ns, double dradius): Vicsek_consensus(ii, vzero, ns) {
    detection_radius2 = dradius * dradius ;
}


int Vicsek_prey::sense_danger(Agent* ag, int num_threats, Agent* threats, double* new_vel) {
    int it , i;
    double disp[DIM] ;
    double dist2 ;
    for(it=0 ; it<num_threats ; it++){
        dist2 = inter->g->distance2(threats[it].get_pos() , ag->get_pos() ) ;
        if(dist2 < detection_radius2){
            inter->g->displacement(threats[it].get_pos() , ag->get_pos(), disp) ;
            for(i=0 ; i<DIM ; i++)
                new_vel[i] = disp[i] * v0 / sqrt(dist2) ;
            return 1 ;
        }
    }
    return 0 ;
}

/*
 * Vicsek Predator
 */
int Vicsek_predator::sense_victims(Agent* pred, int num_agents, Agent* ags){
    int ia , imin;
    double tmp, mindist ;
    imin = 0 ;
    mindist = inter->g->distance2( pred->get_pos(), ags[0].get_pos() ) ;
    for(ia=1; ia<num_agents; ia++){
        tmp = inter->g->distance2( pred->get_pos(), ags[ia].get_pos() ) ;
        if(tmp < mindist){
            mindist = tmp ;
            imin = ia ;
        }
    }
    return imin ;
}

int Vicsek_predator::hunt(Agent* pred, Agent* prey, double dt){
    double disp[DIM] ;
    double dist2 ;
    int i ;
    double* pos = pred->get_pos() ;
    double* vel = pred->get_vel() ;

    inter->g->displacement( pos, prey->get_pos() , disp ) ;
    dist2 = inter->g->length2(disp) ;
    for(i=0; i<DIM; i++)
        vel[i] = disp[i] * v0 / sqrt(dist2) ;

    if(dist2 < v0 * v0 * dt * dt ){
        /* deep copy position, dont copy the pointer! */
        for(i=0; i<DIM; i++)
            pos[i] = prey->get_pos()[i] ;
        return 1 ;
    }else{
        pred->move(dt) ;
        return 0 ;
    }
}
