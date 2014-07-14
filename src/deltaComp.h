/***** deltaComp.h ********************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Fri Jul 11 09:49:15 2014
 **************************************************/
#ifndef DELTACOMP
#define DELTACOMP

#include "profileTree.h"

/* define container for parameters */
typedef struct deltaParam{
  char type;   /* parameter type */
  double de;
  double dLo;
  double dUp;
  double l;
  ProfilePairs *pp;
} DeltaParam;

#endif
