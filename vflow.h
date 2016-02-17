#ifndef VELOCITY
#define VELOCITY

#include <ctime>
#include "vectorXY.h"

int loadvgrid(struct tm rdate,int vfield);
int loadvflow(struct tm seeddate, int ntime, int vfield);

int getvflow(double t,vectorXY point, vectorXY *vint, int vfield);

#endif
