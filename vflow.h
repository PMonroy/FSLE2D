#ifndef VELOCITY
#define VELOCITY

#include <ctime>
#include "vectorXY.h"

/* Global Variables*/
extern int nvlon, nvlat;
extern double *vlon, *vlat;

int loadvgrid(struct tm rdate,int vfield);
void freevgrid(void);

int loadvflow(struct tm seeddate, int ntime, int vfield);
void freevflow(int ntime);

int getvflow(double t,vectorXY point, vectorXY *vint, int vfield);

#endif
