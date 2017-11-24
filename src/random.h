/* Interface to generate random numbers in libspp.
 * spp_set_seed must be called before getting random
 * numbers, or one will get a segmentation fault
 * due to an unallocated pointer.
 *
 * The implementation of random numbers is taken from
 *      http://people.sc.fsu.edu/~jburkardt/c_src/ziggurat/ziggurat.html
 *  (modifications: remove several functions and change float->double)
 *
 * Reference:
 *      George Marsaglia, Wai Wan Tsang,
 *      The Ziggurat Method for Generating Random Variables,
 *      Journal of Statistical Software,
 *      Volume 5, Number 8, October 2000, seven pages.
 */



/* Set the seed value and execute
 * any pre-computation needed to set
 * up the random number generator.
 * This function must be called before
 * using any of the spp_random_*.
 */
void spp_set_seed(long int s) ;
/* Return a normal (gaussian) random
 * number with mean 0 and standard deviation 1.
 * Requires the previous execution of spp_set_seed().
 */
double spp_random_normal() ;
/* Return a uniform random number
 * between 0 and 1.
 * Requires the previous execution of spp_set_seed().
 */
double spp_random_uniform() ;
/* Set the components of *vec* to a set of
 * i.i.d. random normal (gaussian) numbers
 * with mean 0 and standard deviation 1.
 */
void spp_random_normal_vector(double* vec) ;
/* Set the array *vec* to a random vector with
 * norm *norm* and homogeneous angular distribution.
 */
void spp_random_vector(double* vec, double norm) ;
