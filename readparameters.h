#ifndef DATE
#define DATE

#include <ctime>
#include "vectorXY.h"

extern int verbose;
extern int vfield;
extern vectorXY domainmin;
extern vectorXY intergrid;
extern vectorXY domainmax;
extern struct tm seeddate;
extern double  intstep;
extern int tau;
extern double deltamax;

int GetcmdlineParameters(int narg,char ** cmdarg, string *fnameparams);
int GetfileParameters(string nfileparameters);

#endif
