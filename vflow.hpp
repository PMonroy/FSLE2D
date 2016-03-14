#ifndef VELOCITY
#define VELOCITY

#include <ctime>
#include "vectorXY.hpp"
#include "vectorIJ.hpp"

int LoadVGrid(struct tm rdate, vectorXY domainmin, vectorXY domainmax, vectorXY meanvel, double duration);
int LoadVFlow(struct tm seeddate, int ntime);
void FreeVFlow(unsigned int ntime);

int GetVFlow(double t,vectorXY point, vectorXY *vint);

int GetVPartialDeriv(double t,vectorXY point, vectorIJ dir, vectorXY *dvdr);

#endif
