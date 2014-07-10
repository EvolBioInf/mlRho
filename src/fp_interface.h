#ifndef INTERFACE
#define INTERFACE

#define DEFAULT_N "profileDb"
#define DEFAULT_C 4
#define DEFAULT_B 1024

/* define argument container */
typedef struct args{
  char h;   /* help message? */
  char v;   /* version message? */
  char e;   /* error message? */
  char r;   /* profiles only? */
  char *n;  /* name of ourput file */
  int b;    /* buffer size */
  int c;    /* minimum coverage */
  char **inputFiles;
  int numInputFiles;
} Args;

Args *getArgs(int argc, char *argv[]);
void printUsage(char *version);
void printSplash(char *version);

#endif
