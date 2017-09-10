#include <math.h>
class Agent ;

/*
 * Abstract Geometry class used as a template
 * for different Geometry implementations.
 * Every implementation defines a geometry
 * where the agents move, and must implement:
 *
 *      displacement: a rule to determine
 *          the displacement (vector) between
 *          two points.
 *      distance2: a rule to determine the
 *          distance^2 (scalar) between two points.
 */

class Geometry {
    public:
        /* pure virtual, must be implemented */
        virtual void  displacement(double* x0, double* x1, double* dis) = 0 ;
        /* pure virtual, must be implemented */
        virtual double distance2(double* x0, double* x1) = 0 ;
        /* Norm^2 of a vector *vect*.
         * Returns sum_i vect[i]*vect[i].
         */
        double length2(double* vect) ;
        /* Size of the computation box.
         * May be ignored by some implementations.
         */
        double L ;
} ;

/*
 * Abstract Interaction class used as a template
 * for different Interaciton implementations.
 * Each implementation defines a way to determine
 * which agents are connected or "neighbors",
 * and must implement:
 *
 *      get_neighbors: select the neighbors of
 *          *a0* out of a list of *n_agents* agents
 *          stored in *ags*. Store a list of pointers
 *          to the neighbors (NOT a copy of them)
 *          in *neis*.
 *          Return the number of neighbors found.
 *      is_neighbor: return 1 if *a1* is a neighbor
 *          of *a0*, 0 otherwise.Note that this relation
 *          is not symmetric in general.
 *          A call to look_around() may be needed
 *          before using this function.
 *
 */

class Interaction {
    public:
        /* pure virtual, must be implemented */
        virtual int get_neighbors(Agent* a0 , int n_agents , Agent* ags, Agent** neis) = 0 ;
        /* pure virtual, must be implemented */
        virtual int is_neighbor(Agent* a0 , Agent* a1) = 0 ;
        /* Auxiliary function designed for interactions that
         * are not local, meaning interactions where
         * is_neighbor(a,b) does NOT only depend on *a*
         * and *b*, such as topologic interaction.
         * This function can be called to setup whatever
         * is needed for the interaction to work as a local one.
         * By default it does nothing.
         */
        virtual void look_around(Agent* a0 , int n_agents , Agent* ags) {};
        /* Geometry used to measure distances between agents.
         */
        Geometry* g ;
} ;

 //Implementations of Geometry

/* Cartesian geometry, the displacement
 * between two points is the vector
 * difference, and the distance2 is the
 * norm^2 of that vector.
 * This geometry ignores the box size.
 */
class Cartesian: public Geometry {
    public:
        Cartesian(double l) ;
        /* Stores the displacement from *x0* to *x1* in *dis*,
         *      dis = x1 - x0
         */
        void  displacement(double* x0, double* x1, double* dis) ;
        /* Returns the square of the distance between *x0* and *x1*.
         *  distance2 = |dis|^2 = sum_i |x1_i - x0_i|^2
         */
        double distance2(double* x0, double* x1) ;
} ;

/* Same as Cartesian but taking into account
 * periodic boundary conditions of the box.
 * The same boundary conditions apply in
 * all the dimensions.
 * The displacement between two points may
 * be shorter by "going through the wall",
 * i.e. using a mirror copy of either of
 * the points.
 */
class CartesianPeriodic: public Geometry {
    public:
        CartesianPeriodic(double l) ;
        /* The displacement between *x0* and *x1* is
         * stored in *dis*. The displacement between *x0*
         * and *x1* is defined as the vector difference
         * between *x0* or any of its mirror copies
         * and *x1* or any of its mirror copies such that
         * it yields the minimal distance2.
         * In practice, it means to compute the x0-x1
         * for each axis and if the number is bigger
         * than L/2 substract L from it, and if it
         * smaller than -L/2 add L to it.
         */
        void  displacement(double* x0, double* x1, double* dis) ;
        /* Norm2 of displacement.
         */
        double distance2(double* x0, double* x1) ;
} ;

// Implementations of Interaction

/*
 *  Metric interaction: two agents are neighbors
 *  if the distance2 between their positions
 *  (defined by the geometry) is less than or
 *  equal to *rad2*.
 *  For clarity, a radius *r* is given at
 *  initialization which is then squared,
 *  i.e. rad2 = r*r .
 *  This interaction is local and symmetric.
 *
 */
class Metric : public Interaction {
    public:
        Metric() {};
        /* Initialization reads a radius *r*
         * and stores its square in *rad2*.
         */
        Metric(double r , Geometry* g) ;
        /* Copy A POINTER to all the
         * agents in *ags* whose distance2
         * to *a0* is less than or equal to *rad2*
         * into *neis*.
         * Return the number of neighbors found.
         */
        int get_neighbors(Agent* a0 , int n_agents , Agent* ags, Agent** neis) ;
        /* Return 1 if the distance2 between *a0*
         * and *a1* is less than or equal to *rad2*.
         */
        int is_neighbor(Agent* a0 , Agent* a1) ;
        /* Return the interaction radius. */
        double radius(){return sqrt(rad2);} ;
    private:
        double rad2 ;
} ;

/*
 * Topologic interaction: an agent's neighbors
 * are the *k*-closest agents as determined
 * by distance2 from the geometry.
 * This interaction is NOT local as in NOT
 * symmetric. Because it is not local, a
 * call to look_around() is required before
 * using is_neighbor(a0,a1) for each a0.
 *
 * This class requires memory to store the
 * distances to the agents in *dists2*.
 * This memory is NOT allocated on initialization.
 * Also requires a quickselect implementation
 * to sort the distances and return the k-th
 * smallest.
 *
 */
class Topologic : public Interaction {
    public:
        Topologic() {};
        Topologic(int k , Geometry* g, double* dd) ;
        /* Copy A POINTER to the *k*
         * agents in *ags* closer to *a0* into *neis*.
         * Return the number of neighbors found
         * (should always be min(k,n_agents) ).
         */
        int get_neighbors(Agent* a0 , int n_agents , Agent* ags, Agent** neis) ;
        /* Return 1 if *a1* is one of the *k*
         * closest agents to *a0*.
         * WARNING: This requires a call to
         * look_around() for each *a0*.
         */
        int is_neighbor(Agent* a0 , Agent* a1) ;
        /* Computes the effective radius *rad2*
         * of the agent *a0* (the distance2
         * to the k-th closest agent) so that *is_neighbor*
         * can work as in the Metric case.
         */
        void look_around(Agent* a0 , int n_agents , Agent* ags) ;
        /* Return the outdegree = the number of neighbors. */
        int outdegree(){return k;} ;
    private:
        int k ;
        double rad2 ;
        double* dists2 ;
} ;

/*
 * NoInteraction interaction: nobody
 * interacts with anybody. Similar behavior
 * can be achieved with a Metric of r=0
 * of Topologic of k=0 but this is much
 * faster.
 * Beware of incosistency, a0 is not neighbor
 * of itself according to is_neighbor()
 * but it is returned in get_neighbors().
 * This is a small price to pay for simplicity
 * and speed: get_neighbors() needs to
 * return the agent itself for him to
 * maintain current speed in the Vicsek
 * protocol, but adding a if to is_neighbor
 * to check if a0==a1 adds computation
 * cost for no real benefit.
 */
class NoInteraction : public Interaction {
    public:
        NoInteraction( Geometry* g) ;
        int get_neighbors(Agent* a0 , int n_agents , Agent* ags, Agent** neis) ;
        int is_neighbor(Agent* a0 , Agent* a1) ;
} ;

// Other functions

/* Return the *k*-th smaller value
 * out of the *n* values stored in *arr*.
 * Used by Topologic::look_around().
 * Based on the algorithm described
 * in "Numerical Recipes in C".
 *
 * The array *arr* IS modified by
 * this function.
 */
double quickselect(double *arr, int n, int k) ;
