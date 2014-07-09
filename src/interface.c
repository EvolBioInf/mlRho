/***** interface.c ************************************************
 * Description: Routine for gathering arguments from the command
 *              line.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Sun Jun 20 13:12:10 2004.
 *****************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <unistd.h>
#include "interface.h"
#include "eprintf.h"

Args *args;

Args *getArgs(int argc, char *argv[]){
  int c;
  char *optString = "P:E:D:t:s:i:hvM:lLm:S:n:Iof:c";

  args = (Args *)emalloc(sizeof(Args));
  args->P = INI_PI;
  args->E = INI_EPSILON;
  args->D = INI_DELTA;
  args->t = THRESHOLD;
  args->n = DEFAULT_N;
  args->s = STEP_SIZE;
  args->S = DEFAULT_S;
  args->i = MAX_IT;
  args->o = 0;
  args->I = 0;
  args->M = INT_MAX;
  args->c = 0;
  args->m = 1;
  args->l = 0;
  args->L = 0;
  args->h = 0;
  args->e = 0;
  args->v = 0;

  c = getopt(argc, argv, optString);
  while(c != -1){
    switch(c){
    case 'P':                           /* initial value of pi */
      args->P = atof(optarg);
      break;
    case 'E':                           /* initial error rate, epsilon */
      args->E = atof(optarg);
      break;
    case 'c':                           /* corrected diversity computation */
      args->c = 1;
      break;
    case 'D':                           /* initial disequilibrium coefficient, delta */
      args->D = atof(optarg);
      break;
    case 'n':                           /* name of database */
      args->n = optarg;
      break;
    case 'f':                           /* fraction of likelihood weight included in disequilibrium analysis */
      args->f = atof(optarg);
      break;
    case 'o':                           /* High memory mode */
      args->o = 1;
      break;
    case 't':                           /* simplex size threshold */
      args->t = atof(optarg);
      break;
    case 's':                           /* step size in ML computation */
      args->s = atof(optarg);
      break;
    case 'S':                           /* step size in LD computation */
      args->S = atoi(optarg);
      break;
    case 'i':                           /* maximum number of iterations */
      args->i = atoi(optarg);
      break;
    case 'I':                           /* print likelihood values to file */
      args->I = 1;
      break;
    case 'M':                           /* maximum distance investigated in disequilibrium analysis */
      args->M = atoi(optarg);
      break;
    case 'm':                           /* minimum distance investigated in disequilibrium analysis */
      args->m = atoi(optarg);
      break;
    case 'l':                           /* estimate delta */
      args->l = 1;
      break;
    case 'L':                           /* lump "step" distance classes? */
      args->L = 1;
      break;
    case 'v':                           /* print program information */
      args->v = 1;
      break;
    case '?':                           /* fall-through is intentional */
    case 'h':                           /* print help */
      args->h = 1;
      break;
    default:
      printf("# unknown argument: %c\n",c);
      args->e = 1;
      return args;
    }
    c = getopt(argc, argv, optString);
  }
  return args;
}


void printUsage(){
  printf("Usage: mlRho [options] [inputFile(s)]\n");
  printf("Maximum likelihood estimation of population mutation, recombination, and disequilibrium measures\n");
  printf("Example: mlRho -n test -M 0\n");
  printf("standard options:\n");
  printf("\t[-n <FILE> name of database created using formatPro; default: %s]\n",DEFAULT_N);
  printf("\t[-m <NUM> minimum distance analyzed in rho computation; default: 1]\n");
  printf("\t[-M <NUM> maximum distance analyzed in rho computation; default: all]\n");
  printf("\t[-S <NUM> step size in rho computation; default: %d]\n",DEFAULT_S);
  printf("\t[-f <NUM> fraction of likelihood weight included in LD analysis; default: %.2f]\n",DEFAULT_F);
  printf("\t[-o high memory mode; may be faster for disequilibrium analysis]\n");
  printf("\t[-I write likelihoods to file; default: likelihoods not written to file]\n");
  printf("\t[-L lump -S distance classes; default: no lumping]\n");
  printf("\t[-c corrected diversity measure according to Lynch (2008), p. 2412; default: uncorrected]\n");
  printf("\t[-v print program version, etc. and exit]\n");			     
  printf("\t[-h print this help message and exit]\n");
  printf("extra options:\n");
  printf("\t[-P <NUM> initial theta value; default: %10.3e]\n",INI_PI);
  printf("\t[-E <NUM> initial epsilon value; default: %10.3e]\n",INI_EPSILON);
  printf("\t[-D <NUM> initial delta value; default: %10.3e]\n",INI_DELTA);
  printf("\t[-t <NUM> simplex size threshold; default: %10.3e]\n",THRESHOLD);
  printf("\t[-s <NUM> size of first step in ML estimation; default: %10.3e]\n",STEP_SIZE);
  exit(0);
}

void printSplash(){
  const char *str = {
    "***********************************************************\n"
    "*                  mlRho, version " VERSION "                     *\n"
    "*    ML Estimation of Mutation and Recombination Rate     *\n"
    "*  Bernhard Haubold, Peter Pfaffelhuber, Michael Lynch    *\n"
    "*---------------------------------------------------------*\n"
    "* REFERENCES                                              *\n"
    "* 1) Lynch, M. (2008). Estimation of nucleotide           *\n"
    "*    diversity, disequilibrium coefficients, and mutation *\n"
    "*    rates from high-coverage genome-sequencing projects. *\n"
    "*    Mol. Biol. Evol., 25:2409--2419.                     *\n"
    "* 2) Haubold, B., Pfaffelhuber, P. and Lynch, M.          *\n"
    "*    (2009). mlRho - A program for estimating the pop-    *\n"
    "*    ulation mutation and recombination rates from shot-  *\n"
    "*    gun-sequenced diploid genomes. Mol. Ecol., 19:277-   *\n"
    "*    284.                                                 *\n"
    "* 3) Thota, A., Haubold, B., Michael, S, Doak,T.,         *\n"
    "*    Xu, S., and Henschel, R. (2013). Making campus       *\n"
    "*    bridging work for researchers: a case study with     *\n"
    "*    mlRho. In: Proceedings of the Conference on Extreme  *\n"
    "*    Science and Engineering Discovery Environment:       *\n"
    "*    Gateway to Discovery, 17:1-17:8. ACM, New York.      *\n"
    "*                                                         *\n"  
    "* CONTACT                                                 *\n"
    "* Code maintained by Bernhard Haubold,                    *\n"
    "* haubold@evolbio.mpg.de                                  *\n"
    "*                                                         *\n"
    "* LICENSE                                                 *\n"
    "* This software is distributed under the GNU General      *\n"
    "* Public License. You should have received a copy         *\n"
    "* of the licence together with this software. If          *\n"
    "* not, see http://www.gnu.org/licenses/.                  *\n"
    "***********************************************************\n"
  };

  printf("%s", str);
  exit(0);
}

     
