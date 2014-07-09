/***** filterPro.c ********************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Sun Jul  7 17:18:48 2013
 **************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "profile.h"
#include "mlComp.h"
#include "eprintf.h"

int ascCompDouble(const void *a, const void *b);
int descCompDouble(const void *a, const void *b);

char *filterPro(double pi, double inclFract){
	Profile *profiles;
	size_t i, numProfiles;
	double *lik, *originalLik, piComp, sum;
	double minLik, minSum;
	char *incl;

	numProfiles = (size_t) getNumProfiles();

	incl = (char *) emalloc(numProfiles*sizeof(char));
	if(inclFract == 1.0){
		memset( incl, 1, numProfiles);
		return incl;
	}

	profiles = getProfiles();

	double const *lOnes = getLones();
	double const *lTwos = getLtwos();

	lik = (double *) emalloc(numProfiles*sizeof(double));
	originalLik = (double *) emalloc(numProfiles*sizeof(double));

	piComp = 1.0 - pi;
	sum = 0.0;
	for(i = 0; i < numProfiles; i++){
		lik[i] = -log(lOnes[i]*piComp + lTwos[i]*pi) * profiles[i].n;
		originalLik[i] = lik[i];
		sum += lik[i];
	}

	qsort( lik, numProfiles, sizeof(double), descCompDouble);

	minSum = sum * inclFract;
	sum = 0.0;
	for(i = 0; sum < minSum; i++){
		sum += lik[i];
	}

	minLik = lik[i-1];
	for(i = 0; i < numProfiles; i++){
		incl[i] = originalLik[i] >= minLik ? 1 : 0;
	}

	free(originalLik);
	free(lik);
	return incl;
}

// This function is for ascending sorting. The original implementaion
// of compDouble was wrong.
int ascCompDouble(const void *a, const void *b){
	const double aa = *((double *) a);
	const double bb = *((double *) b);

	if( aa > bb ) return 1;
	if( aa < bb ) return -1;
	return 0;
}

int descCompDouble(const void *a, const void *b){
	return -ascCompDouble(a,b);
}
