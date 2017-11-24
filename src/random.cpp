#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include "random.h"

double r4_nor ( uint32_t *jsr, uint32_t kn[128], double fn[128], double wn[128] );
void r4_nor_setup ( uint32_t kn[128], double fn[128], double wn[128] );
double r4_uni ( uint32_t *jsr );
uint32_t shr3_seeded ( uint32_t *jsr );

/* The seed and the pointer to it are
 * defined separately so that the pointer
 * can be initialized to NULL. This is
 * done to make sure that trying to
 * generate random numbers before calling
 * spp_set_seed (and doing the required
 * setup e.g. r4_nor_setup) will fail.
 */

uint32_t* spp_seed_ptr = NULL; // must call set_seed to use.
uint32_t spp_seed ;
uint32_t kn[128] ;
double fn[128] ;
double wn[128] ;

void spp_set_seed(long int s){
    r4_nor_setup(kn, fn, wn);
    spp_seed = (uint32_t) s ;
    spp_seed_ptr = &spp_seed ;
}

double spp_random_normal(){
    return r4_nor(spp_seed_ptr, kn, fn, wn) ;
}

double spp_random_uniform(){
    return r4_uni(spp_seed_ptr) ;
}

void spp_random_normal_vector(double* vec){
    /* Note that a vector whose components are
     * normally distributed has a normally distributed
     * norm and homogeneous angular distribution.
     * This is NOT true in general (i.e. a vector
     * with uniformly distributed components does
     * not have uniformly distributed norm.)
     */
    int i ;
    for(i=0; i<DIM; i++)
        vec[i] = r4_nor(spp_seed_ptr, kn, fn, wn) ;
}

void spp_random_vector(double* vec, double norm){
    /* To fix the norm, first generate a normally
     * distributed vector and then re-scale it.
     */
    int i ;
    double v2 = 0.0 ;
    for(i=0; i<DIM; i++){
        vec[i] = r4_nor(spp_seed_ptr, kn, fn, wn) ;
        v2 += vec[i]*vec[i] ;
    }
    for(i=0; i<DIM; i++)
        vec[i] *= norm/sqrt(v2) ;
}


/*
 * The functions below this are taken from
 *      http://people.sc.fsu.edu/~jburkardt/c_src/ziggurat/ziggurat.html
 *  (modifications: remove several functions and change float->double)
 *
 * Reference:
 *      George Marsaglia, Wai Wan Tsang,
 *      The Ziggurat Method for Generating Random Variables,
 *      Journal of Statistical Software,
 *      Volume 5, Number 8, October 2000, seven pages.
 *
 */


double r4_nor ( uint32_t *jsr, uint32_t kn[128], double fn[128], double wn[128] )

/******************************************************************************/
/*
  Purpose:

    R4_NOR returns a normally distributed single precision real value.

  Discussion:

    The value returned is generated from a distribution with mean 0 and
    variance 1.

    The underlying algorithm is the ziggurat method.

    Before the first call to this function, the user must call R4_NOR_SETUP
    to determine the values of KN, FN and WN.

    Thanks to Chad Wagner, 21 July 2014, for noticing a bug of the form
      if ( x * x <= y * y );   <-- Stray semicolon!
      {
        break;
      }

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 July 2014

  Author:

    John Burkardt

  Reference:

    George Marsaglia, Wai Wan Tsang,
    The Ziggurat Method for Generating Random Variables,
    Journal of Statistical Software,
    Volume 5, Number 8, October 2000, seven pages.

  Parameters:

    Input/output, uint32_t *JSR, the seed.

    Input, uint32_t KN[128], data computed by R4_NOR_SETUP.

    Input, double FN[128], WN[128], data computed by R4_NOR_SETUP.

    Output, double R4_NOR, a normally distributed random value.
*/
{
  int hz;
  uint32_t iz;
  const double r = 3.442620;
  double value;
  double x;
  double y;

  hz = ( int ) shr3_seeded ( jsr );
  iz = ( hz & 127 );

  if ( fabs ( hz ) < kn[iz] )
  {
    value = ( double ) ( hz ) * wn[iz];
  }
  else
  {
    for ( ; ; )
    {
      if ( iz == 0 )
      {
        for ( ; ; )
        {
          x = - 0.2904764 * log ( r4_uni ( jsr ) );
          y = - log ( r4_uni ( jsr ) );
          if ( x * x <= y + y )
          {
            break;
          }
        }

        if ( hz <= 0 )
        {
          value = - r - x;
        }
        else
        {
          value = + r + x;
        }
        break;
      }

      x = ( double ) ( hz ) * wn[iz];

      if ( fn[iz] + r4_uni ( jsr ) * ( fn[iz-1] - fn[iz] )
        < exp ( - 0.5 * x * x ) )
      {
        value = x;
        break;
      }

      hz = ( int ) shr3_seeded ( jsr );
      iz = ( hz & 127 );

      if ( fabs ( hz ) < kn[iz] )
      {
        value = ( double ) ( hz ) * wn[iz];
        break;
      }
    }
  }

  return value;
}
/******************************************************************************/

void r4_nor_setup ( uint32_t kn[128], double fn[128], double wn[128] )

/******************************************************************************/
/*
  Purpose:

    R4_NOR_SETUP sets data needed by R4_NOR.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 October 2013

  Author:

    John Burkardt

  Reference:

    George Marsaglia, Wai Wan Tsang,
    The Ziggurat Method for Generating Random Variables,
    Journal of Statistical Software,
    Volume 5, Number 8, October 2000, seven pages.

  Parameters:

    Output, uint32_t KN[128], data needed by R4_NOR.

    Output, double FN[128], WN[128], data needed by R4_NOR.
*/
{
  double dn = 3.442619855899;
  int i;
  const double m1 = 2147483648.0;
  double q;
  double tn = 3.442619855899;
  const double vn = 9.91256303526217E-03;

  q = vn / exp ( - 0.5 * dn * dn );

  kn[0] = ( uint32_t ) ( ( dn / q ) * m1 );
  kn[1] = 0;

  wn[0] = ( double ) ( q / m1 );
  wn[127] = ( double ) ( dn / m1 );

  fn[0] = 1.0;
  fn[127] = ( double ) ( exp ( - 0.5 * dn * dn ) );

  for ( i = 126; 1 <= i; i-- )
  {
    dn = sqrt ( - 2.0 * log ( vn / dn + exp ( - 0.5 * dn * dn ) ) );
    kn[i+1] = ( uint32_t ) ( ( dn / tn ) * m1 );
    tn = dn;
    fn[i] = ( double ) ( exp ( - 0.5 * dn * dn ) );
    wn[i] = ( double ) ( dn / m1 );
  }

  return;
}
/******************************************************************************/

double r4_uni ( uint32_t *jsr )

/******************************************************************************/
/*
  Purpose:

    R4_UNI returns a uniformly distributed real value.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 October 2013

  Author:

    John Burkardt

  Reference:

    George Marsaglia, Wai Wan Tsang,
    The Ziggurat Method for Generating Random Variables,
    Journal of Statistical Software,
    Volume 5, Number 8, October 2000, seven pages.

  Parameters:

    Input/output, uint32_t *JSR, the seed.

    Output, double R4_UNI, a uniformly distributed random value in
    the range [0,1].
*/
{
  uint32_t jsr_input;
  double value;

  jsr_input = *jsr;

  *jsr = ( *jsr ^ ( *jsr <<   13 ) );
  *jsr = ( *jsr ^ ( *jsr >>   17 ) );
  *jsr = ( *jsr ^ ( *jsr <<    5 ) );

  value = fmod ( 0.5
    + ( double ) ( jsr_input + *jsr ) / 65536.0 / 65536.0, 1.0 );

  return value;
}
/******************************************************************************/

uint32_t shr3_seeded ( uint32_t *jsr )

/******************************************************************************/
/*
  Purpose:

    SHR3_SEEDED evaluates the SHR3 generator for integers.

  Discussion:

    Thanks to Dirk Eddelbuettel for pointing out that this code needed to
    use the uint32_t data type in order to execute properly in 64 bit mode,
    03 October 2013.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 October 2013

  Author:

    John Burkardt

  Reference:

    George Marsaglia, Wai Wan Tsang,
    The Ziggurat Method for Generating Random Variables,
    Journal of Statistical Software,
    Volume 5, Number 8, October 2000, seven pages.

  Parameters:

    Input/output, uint32_t *JSR, the seed, which is updated
    on each call.

    Output, uint32_t SHR3_SEEDED, the new value.
*/
{
  uint32_t value;

  value = *jsr;

  *jsr = ( *jsr ^ ( *jsr <<   13 ) );
  *jsr = ( *jsr ^ ( *jsr >>   17 ) );
  *jsr = ( *jsr ^ ( *jsr <<    5 ) );

  value = value + *jsr;

  return value;
}
/******************************************************************************/
