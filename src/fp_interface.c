/***** interface.c ************************************************
 * Description: Routine for gathering arguments from the command
 *              line.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Sun Jun 20 13:12:10 2004.
 *****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "fp_interface.h"
#include "eprintf.h"

Args *args;

Args *getArgs(int argc, char *argv[]){
  char c;
  char *optString = "hvn:rb:c:";

  args = (Args *)emalloc(sizeof(Args));
  args->n = DEFAULT_N;
  args->b = DEFAULT_B;
  args->c = DEFAULT_C;
  args->r = 0;
  args->h = 0;
  args->v = 0;
  args->e = 0;

  c = getopt(argc, argv, optString);
  while(c != -1){
    switch(c){
    case 'r':                           /* profiles only? */
      args->r = 1;
      break;
    case 'b':                           /* buffer size */
      args->b = atoi(optarg);
      break;
    case 'c':                           /* minimum coverage */
      args->c = atoi(optarg);
      break;
    case 'n':                           /* name of output files */
      args->n = optarg;
      break;
    case '?':                           /* fall-through is intentional */
    case 'h':                           /* print help */
      args->h = 1;
      break;
    case 'v':                           /* print version */
      args->v = 1;
      break;
    default:
      printf("# unknown argument: %c\n",c);
      args->e = 1;
      return args;
    }
    c = getopt(argc, argv, optString);
  }
  args->inputFiles = argv + optind;
  args->numInputFiles = argc - optind;
  return args;
}


void printUsage(char *version){
  printf("Usage: %s [options] [inputFiles]\n",progname());
  printf("Format profiles for subsequent analysis with proStats and/or mlRho\n");
  printf("formatPro myProfiles1.pro myProfiles2.pro\n");
  printf("Options:\n");
  printf("\t[-n <FileName>; default: %s]\n",DEFAULT_N);
  printf("\t[-b <NUM> buffer size; default: %d]\n",DEFAULT_B);
  printf("\t[-c <NUM> minimum coverage; default: %d]\n",DEFAULT_C);
  printf("\t[-r create only profiles file; default: profiles, contigs, and positions]\n");
  printf("\t[-h print this help message and exit]\n");
  printf("\t[-v print program information and exit]\n");
  exit(0);
}

void printSplash(char *version){
  printf("%s %s\n",progname(),version);
  printf("Written by Bernhard Haubold.\n");
  printf("Distributed under the GNU General Public License.\n");
  printf("Please send bug reports to haubold@evolbio.mpg.de\n");
  exit(0);
}
