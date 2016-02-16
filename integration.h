
int RK4(double t0, double intstep, vectorXY *point, int (*velocity)(double t,vectorXY point, vectorXY *vint, int vfield), int vfield);
