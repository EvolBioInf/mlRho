/***** inspectPro.h *******************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Wed Nov 14 17:06:14 2012
 **************************************************/
#ifndef INSPECTPRO
#define INSPECTPRO

typedef struct result{
  double pi;   
  double ee;
  double de;   /* delta */
  double rh;   /* rho */
  double rhoFromDelta;
  double l;    /* likelihood */		 
  double pLo;  /* lower bound of pi */
  double pUp;  /* upper bound of pi */
  double eLo;  /* lower bound of epsilon */
  double eUp;  /* upper bound of epsilon */
  double dLo;  /* lower bound of delta */
  double dUp;  /* upper bound of delta */
  double rLo;  /* lower bound of rho */
  double rUp;  /* upper bound of rho */
  int i;       /* number of iterations */
  char type;   /* parameter type */
}Result;

#endif
