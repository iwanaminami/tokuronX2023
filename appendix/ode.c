#include <R.h>
/* a trick to keep up with the parameters */

static double parms[6];

#define lambda parms[0]
#define beta parms[1]
#define d parms[2]
#define p parms[3]
#define delta parms[4]
#define c parms[5]


/* initializers */
void initparms( void (* odeparms)(int *, double *) ) {
  int N = 6;
  odeparms(&N,parms);
}
/* names for states and derivatives */
#define Tp var[0]
#define Ip var[1]
#define Vp var[2]

#define dTpdt vardot[0]
#define dIpdt vardot[1]
#define dVpdt vardot[2]


void derivs( int *neq, double *t, double *var, double *vardot, double *varout, int *ip ) {
  if( ip[0]<1 ) {
    error("nout should be at least 1");
  }
  
  
  dTpdt = lambda - beta*Tp*Vp - d*Tp;
  dIpdt = beta*Tp*Vp - delta*Ip;
  dVpdt = p*Ip - c*Vp;
  
  
}
