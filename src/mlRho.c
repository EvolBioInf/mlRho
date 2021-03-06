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
#include "eprintf.h"
#include "interface.h"
#include "ld.h"
#include "profile.h"
#include "profileTree.h"
#include "mlComp.h"

void runAnalysis(Args *args);
void freeMem(Node **profilePairs, int numProfiles);

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
  Result *r;
  int i;
  int numProfiles;
  char *headerPi, *headerDeltaRho, *outStrPi, *outStrDeltaRho;
  char *inclPro; /* indicates whether or not to include a profile */
  double numPos, c, n;
  Profile *profiles;
  ContigDescr *contigDescr;
  FILE *fp;
  Node **profilePairs;

  headerPi = "d\tn\ttheta\t\t\t\tepsilon\t\t\t\t-log(L)\n";
  headerDeltaRho = "d\tn\ttheta\t\t\t\tepsilon\t\t\t\t-log(L)\t\tdelta\t\t\t\trho\n";
  outStrPi = "%d\t%.0f\t%8.2e<%8.2e<%8.2e\t%8.2e<%8.2e<%8.2e\t%8.2e\n";
  outStrDeltaRho = "%d\t%.0f\t\t\t\t\t\t\t\t\t%8.2e\t%8.2e<%8.2e<%8.2e\t%8.2e<%8.2e<%8.2e\n";
  r = newResult();
  /* heterozygosity analysis */
  readProfiles(args->n);
  numProfiles = getNumProfiles();
  profiles = getProfiles();
  if(args->M == 0 && numProfiles)
    printf("%s", headerPi);
  else
    printf("%s", headerDeltaRho);
  if(numProfiles){
    r = estimatePi(profiles,numProfiles,args,r);
    numPos = piComp_getNumPos(profiles, numProfiles);
    if(args->c){
      /* including the correction from Lynch (2008), p. 2412 */
      n = getCoverage();
      c = n*pow(0.5,n-1);
      c = 1-c;
      r->pLo /= c;
      r->pi /= c;
      r->pUp /= c;
    }
    printf(outStrPi,0,numPos,r->pLo,r->pi,r->pUp,r->eLo,r->ee,r->eUp,r->l);
  }
  fflush(NULL);
  /* linkage analysis */
  fp = iniLdAna(args);
  inclPro = filterPro(r->pi,args->f);
  contigDescr = getContigDescr();
  profilePairs = NULL;
  if(args->o)
    readPositions(fp, contigDescr);
  for(i=args->m;i<=args->M;i+=args->S){
    profilePairs = getProfilePairs(numProfiles, inclPro, contigDescr, fp, args, i);
    r = estimateDelta(profilePairs,numProfiles,args,r,i);
    printf(outStrDeltaRho,i,getNumPos(),r->l,r->dLo,r->de,r->dUp,r->rLo,r->rh,r->rUp);
    fflush(NULL);
  }
  fclose(fp);
  if(args->I){ 
    if(args->c){
      /* save uncorrected values */
      n = getCoverage();
      c = n*pow(0.5,n-1);
      c = 1-c;
      r->pLo *= c;
      r->pi *= c;
      r->pUp *= c;
    }
    writeLik(args->n,r);
  }
  free(r);
  freeMem(profilePairs, numProfiles);
}

void freeMem(Node **profilePairs, int numProfiles){
  int i;
  double *lOnes, *lTwos;
  Profile *profiles;
  ContigDescr *cd;

  if(profilePairs){
    for(i=0;i<numProfiles;i++)
      freeTree(profilePairs[i]);
    free(profilePairs);
  }
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
