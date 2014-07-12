/***** mlRho.c ************************************
 * Description: Maximum-likelihood estimation of 
 *   mutation and recombination rates from re-
 *   sequencing data.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Wed Feb 18 16:35:16 2009.
 **************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <assert.h>
#include <omp.h>
#include "eprintf.h"
#include "interface.h"
#include "ld.h"
#include "profile.h"
#include "profileTree.h"
#include "mlComp.h"

void runAnalysis(Args *args);
void freeMem();

int main(int argc, char *argv[]){
  Args *args;
  char *version;

  version = "2.7";
  setprogname2("mlRho");
  args = getArgs(argc, argv);
  if(args->v)
    printSplash(version);
  if(args->h || args->e)
    printUsage(version);
  runAnalysis(args);
  free(args);
  free(progname());
  return 0;
}

void runAnalysis(Args *args){
  Result *r, *piRes, **results;
  int i, *distArr, numRes;
  int numProfiles;
  char *headerPi, *headerDeltaRho, *outStrPi, *outStrDeltaRho;
  char *inclPro; /* indicates whether or not to include a profile */
  double numPos, c, n;
  Profile *profiles;
  long *numPosArr;
  ContigDescr *contigDescr;
  FILE *fp;
  ProfilePairs *profilePairs;


  headerPi = "d\tn\ttheta\t\t\t\tepsilon\t\t\t\t-log(L)\n";
  headerDeltaRho = "d\tn\ttheta\t\t\t\tepsilon\t\t\t\t-log(L)\t\tdelta\t\t\t\trho\n";
  outStrPi = "%d\t%.0f\t%8.2e<%8.2e<%8.2e\t%8.2e<%8.2e<%8.2e\t%8.2e\n";
  outStrDeltaRho = "%d\t%ld\t\t\t\t\t\t\t\t\t%8.2e\t%8.2e<%8.2e<%8.2e\t%8.2e<%8.2e<%8.2e\n";
  piRes = newResult();
  /* heterozygosity analysis */
  readProfiles(args->n);
  numProfiles = getNumProfiles();
  profiles = getProfiles();
  if(args->M == 0 && numProfiles)
    printf("%s", headerPi);
  else
    printf("%s", headerDeltaRho);
  if(numProfiles){
    piRes = estimatePi(profiles,numProfiles,args,piRes);
    numPos = piComp_getNumPos(profiles, numProfiles);
    if(args->c){
      /* including the correction from Lynch (2008), p. 2412 */
      n = getCoverage();
      c = n*pow(0.5,n-1);
      c = 1-c;
      piRes->pLo /= c;
      piRes->pi /= c;
      piRes->pUp /= c;
    }
    printf(outStrPi,0,numPos,piRes->pLo,piRes->pi,piRes->pUp,piRes->eLo,piRes->ee,piRes->eUp,piRes->l);
  }
  fflush(NULL);
  /* linkage analysis */
  fp = iniLdAna(args);
  inclPro = filterPro(piRes->pi,args->f);
  contigDescr = getContigDescr();
  profilePairs = NULL;
  if(args->o)
    readPositions(fp, contigDescr);
  numRes = (args->M - args->m + 1)/args->S + 1;
  results = (Result **)emalloc(numRes*sizeof(Result *));
  numPosArr = (long *)emalloc(numRes*sizeof(long));
  distArr = (int *)emalloc(numRes*sizeof(int));
  r = copyResult(piRes);
  for(i=0;i<numRes;i++)
    results[i] = copyResult(r);
  free(r);
  numRes = 0;
  printf("args->T: %d\n",args->T);
#pragma omp parallel for num_threads(args->T)
  for(i=args->m;i<=args->M;i+=args->S){
    profilePairs = getProfilePairs(numProfiles, inclPro, contigDescr, fp, args, i);
    profilePairs->dist = i;
    r = results[numRes];
    r = estimateDelta(profilePairs,args,r);
    numPosArr[numRes] = profilePairs->numPos;
    distArr[numRes] = i;
    numRes++;
    freeProfilePairs(profilePairs);
  }
  for(i=0;i<numRes;i++){
    r = results[i];
    printf(outStrDeltaRho,distArr[i],numPosArr[i],r->l,r->dLo,r->de,r->dUp,r->rLo,r->rh,r->rUp);
  }
  fclose(fp);
  if(args->I){ 
    if(args->c){
      /* save uncorrected values */
      n = getCoverage();
      c = n*pow(0.5,n-1);
      c = 1-c;
      piRes->pLo *= c;
      piRes->pi *= c;
      piRes->pUp *= c;
    }
    writeLik(args->n,piRes);
  }
  for(i=0;i<numRes;i++)
    free(results[i]);
  free(results);
  free(inclPro);
  free(distArr);
  free(numPosArr);
  free(piRes);
  freeMem();
}

void freeMem(){
  int i;
  double *lOnes, *lTwos;
  Profile *profiles;
  ContigDescr *cd;

  cd = getContigDescr();
  if(cd){
    if(cd->pos){
      for(i=0;i<cd->n;i++)
	free(cd->pos[i]);
      free(cd->pos);
    }
    free(cd->posBuf);
    free(cd->len);
    free(cd);
  }
  lOnes = getLones();
  lTwos = getLtwos();
  profiles = getProfiles();
  free(lOnes);
  free(lTwos);
  free(profiles);
  freeMlComp();
}
