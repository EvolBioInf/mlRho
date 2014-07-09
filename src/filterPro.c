/***** filterPro.c ********************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Sun Jul  7 17:18:48 2013
 **************************************************/
#include <stdio.h>
#include <stdlib.h>
#include "profile.h"
#include "mlComp.h"
#include "eprintf.h"

int compDouble(const void *a, const void *b);

char *filterPro(double pi, double inclFract){
  Profile *profiles;
  int i, numProfiles;
  double *lOnes, *lTwos;
  double *lik, *originalLik, piComp, sum;
  double minLik, minSum;
  char *incl;

  numProfiles = getNumProfiles();

  incl = (char *)emalloc(numProfiles*sizeof(char));
  if(inclFract == 1.0){
    for(i=0;i<numProfiles;i++)
      incl[i] = 1;
    return incl;
  }
  profiles = getProfiles();
  lOnes = getLones();
  lTwos = getLtwos();

  lik = (double *)emalloc(numProfiles*sizeof(double));
  originalLik = (double *)emalloc(numProfiles*sizeof(double));
  piComp = 1. - pi;
  sum = 0.0;
  for(i=0;i<numProfiles;i++){
    lik[i] = -log(lOnes[i]*piComp + lTwos[i]*pi) * profiles[i].n;
    originalLik[i] = lik[i];
    sum += lik[i];
  }
  qsort(lik,numProfiles,sizeof(double),compDouble);
  minSum = sum * inclFract;
  sum = 0;
  i = 0;
  while(sum < minSum)
    sum += lik[i++];
  minLik = lik[i-1];
  for(i=0;i<numProfiles;i++)
    if(originalLik[i] >= minLik)
      incl[i] = 1;
    else
      incl[i] = 0;
  return incl;
}

int compDouble(const void *a, const void *b){
  double x;

  x = *((double *) b) - *((double *) a); 

  if(x > 0)
    return 1;
  else if(x < 0)
    return -1;
  else
    return 0;
}
