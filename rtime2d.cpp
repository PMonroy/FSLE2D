/* fsle2d.cpp*/
#include <cstdio>
#include <cstdlib>
#include <iostream> // cout, endl...  
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>
#include <iomanip> // setfill, setw
#include <string>

using namespace std;

#include "rparameters.hpp"
#include "vectorIJ.hpp"
#include "gridconstruction.hpp" 
#include "constants.hpp"
#include "vflow.hpp" 
#include "integration.hpp"

// Maybe one day I use this:
//int (*velocity)(double ,vectorXY , vectorXY* ); 

string numprintf(int ndigits, int ndecimals, double number);
struct rtime2dParameters {

  const vectorXY domainmin;
  const vectorXY domainmax;
  const vectorXY intergrid;
  const struct tm seeddate;
  const double intstep;
  const int tau;

  // Define a constructor that will load stuff from a configuration file.
  rtime2dParameters(const string & rtime2dParamsFileName)
  :domainmin(getVectorXYParam(rtime2dParamsFileName, "domainmin"))
  ,domainmax(getVectorXYParam(rtime2dParamsFileName, "domainmax"))
  ,intergrid(getVectorXYParam(rtime2dParamsFileName, "intergrid")) 
  ,seeddate(getDateParam(rtime2dParamsFileName, "seeddate")) 
  ,intstep(getDoubleParam(rtime2dParamsFileName, "intstep")) 
  ,tau(getIntParam(rtime2dParamsFileName, "tau"))
{}
};

int main(int argc, char **argv)
{

  /********************************************************************************
   * READ PARAMETERS FROM FILE
   ********************************************************************************/

  string fnameparams;// File name that stores the parameters
  int namefileflag;
  if(GetcmdlineParameters(argc, argv, &fnameparams, &namefileflag)){//Get command line parameters
    return 1;
  }

  const rtime2dParameters rtime2dParams(fnameparams);

  /* VERBOSE */
#ifdef DEBUG
    cout << "FSLE2D parameters from file: ";
    cout << fnameparams <<endl; 
    cout << " domainmin = "<< rtime2dParams.domainmin<<endl;
    cout << " intergrid = " <<rtime2dParams.intergrid<<endl;
    cout << " domainmax = "<< rtime2dParams.domainmax<<endl;
    cout << " seeddate = "<< rtime2dParams.seeddate.tm_mday;
    cout << " " << rtime2dParams.seeddate.tm_mon+1;
    cout << " " << rtime2dParams.seeddate.tm_year<<endl;
    cout << " intstep = "<<rtime2dParams.intstep<<endl ;
    cout << " tau = "<<rtime2dParams.tau<<endl ;
    cout << endl;
#endif

  /********************************************************************************
   * GRID CONSTRUCTION
   ********************************************************************************/

  vector<vectorXY> grid;
  vectorIJ dimension;

  if(MakeRegularGrid(&grid, &dimension, rtime2dParams.domainmin, rtime2dParams.intergrid, rtime2dParams.domainmax)) {
    cout << "Error in contruction grid" << endl;
    return 1;
  }
  unsigned int numgridpoints=grid.size(); 
  vector<int> qflag(numgridpoints, 1);

#ifdef DEBUG
  cout << "GRID CONSTRUCTION: "<<endl; 
  cout << " num. nodes = "<< grid.size()<<endl;
  cout << " dim.i = "<<dimension.i<<endl;
  cout << " dim.j = "<<dimension.j<<endl;
  cout << endl;
#endif

  /********************************************************************************
   * SETUP TIME PARAMETERS
   ********************************************************************************/

  struct tm *inidate = {0};

  double tstart;
  double tend;
  double h;

  int ntime = abs(rtime2dParams.tau);
  int ascnd = rtime2dParams.tau > 0;
  struct tm datebuff = rtime2dParams.seeddate; // Buffer var to avoid error "cast const tm* to tm*"
  time_t seedtime = mktime(&datebuff); // Convert date to time in seconds (UTC) 

  if(ascnd){
    time_t initime = seedtime;
    inidate = gmtime(&initime);
    tstart = 0.0;
    tend = (double) ntime;
    h= rtime2dParams.intstep;
  }
  else{
    time_t initime = seedtime - ntime*secondsday;
    inidate = gmtime(&initime); 
    tend = 0.0;
    tstart = (double) ntime;
    h=-1.0*rtime2dParams.intstep;
  }

  /********************************************************************************
   * LOAD VELOCITY FIELD
   ********************************************************************************/

#ifdef DEBUG
  cout << "Loading velocity grid:" << endl;
#endif

 vectorXY meanvel(0.2,0.2);
 if(LoadVGrid(rtime2dParams.seeddate,
	      rtime2dParams.domainmin,
	      rtime2dParams.domainmax, 
	      meanvel, 
	      abs(rtime2dParams.tau))!=0){//Load velocity grid
   cout << "Error reading velocity grid" << endl;
   return 1;
 }

#ifdef DEBUG  
 cout << "[Complete]" << endl;
 cout << "Loading velocities:" << endl;
#endif

 
 if(LoadVFlow(*inidate,ntime+2)!=0) {  // Load velocity field
   cout << "Error reading velocities"<< endl;
   return 1;
 }
  
#ifdef DEBUG 
  cout << "[Complete]"<<endl;
#endif

  /********************************************************************************
   * TIME LOOP
   ********************************************************************************/
  
  vector<vectorXY > tracer = grid;  
  vector<vectorXY > dvdx;
  vector<vectorXY > dvdy;

  vector<double > ow;
  vector<double > rtime;
  
  vectorIJ dirx(1,0),diry(0,1);

  dvdx.resize(numgridpoints);
  dvdy.resize(numgridpoints);
  ow.resize(numgridpoints);
  rtime.resize(numgridpoints);

  unsigned int step;
  double t;
  unsigned int q;

#ifdef DEBUG//Verbose: Response and exit time calc. 
  cout << "Calculation of retention time:" << endl;
#endif

  for(t=tstart, step=0; ascnd==1?(t<tend):(t>=tend); t+=h,step++){

#ifdef DEBUG//Verbose: show the current step 
    cout << "step=" << step << "(" <<(tend-tstart)/h <<")"<<endl;
#endif

    /* Compute the position of tracer in time t=t+h*/
    for (q=0; q<numgridpoints; q++){
      if(qflag[q]==0 || qflag[q]==-1){//Skip the rest of the loop if q is a non-integrable point 
	  continue;
	}
	
      if(RK4(t, h, &tracer[q], GetVFlow)==1){
	qflag[q]=-1; // q is a non-integrable grid point
      }
    }

    /* Compute the OKUBO-WEISS */
    for (q=0; q<numgridpoints; q++){

      if(qflag[q]==0 || qflag[q]==-1){//Skip the rest of the loop if q is a non-integrable point 
	  continue;
	}
	
      if(GetVPartialDeriv(t,tracer[q], dirx, &dvdx[q])==1){
	qflag[q]=-1;
	continue;
      }
      if(GetVPartialDeriv(t,tracer[q], diry, &dvdy[q])==1){
	qflag[q]=-1;
	continue;
      }
      
      ow[q]=(dvdx[q].x-dvdy[q].y)*(dvdx[q].x -dvdy[q].y)+4.0*(dvdy[q].x*dvdx[q].y);//Okubo-Weiss Parameters
      if(ow[q]>=0.0){
	rtime[q]=step*abs(h);
	qflag[q]=0;
      }
    }
  }

  //free vflow
  FreeVFlow((unsigned int)(ntime+2));

#ifdef DEBUG//Verbose: Success response and exit time calculation
    cout << "[Complete]" << endl;
#endif
  
  /********************************************************************************
   * SAVE RESULTS IN FILES
   ********************************************************************************/

  // Save rtime 2d grid in a file
  string nfilegridrtime2d =     
    "lon"    + numprintf(4,0,rtime2dParams.domainmin.x) 
    + numprintf(4,0,rtime2dParams.domainmax.x)
    + numprintf(4,3,rtime2dParams.intergrid.x)
    + "_lat"       + numprintf(4,0,rtime2dParams.domainmin.y)
    + numprintf(4,0,rtime2dParams.domainmax.y)
    + numprintf(4,3,rtime2dParams.intergrid.y)
    + ".grid";

#ifdef DEBUG
  cout << "Save grid in file: " << nfilegridrtime2d <<endl;
#endif

  if(!ifstream(nfilegridrtime2d.c_str())){ // Check if file exists
    ofstream ofilegridrtime2d(nfilegridrtime2d.c_str());
    for(q=0; q<numgridpoints; q++) {
      ofilegridrtime2d<<grid[q].x<<" "<<grid[q].y<<endl;
    }
    ofilegridrtime2d.close();
  }
  
#ifdef DEBUG
  cout << "[Complete]" << endl;
#endif

  // save fsle 2d values in a file
  string rawname;

  if(namefileflag==1){
    size_t lastdot = fnameparams.find_last_of(".");
    if(lastdot == string::npos){
      rawname="rtime2d_"+fnameparams;
    }
    else{
      rawname="rtime2d_"+fnameparams.substr(0,lastdot);
    }
  }
  else{
    rawname = "rtime2d_lon"    + numprintf(4,0,rtime2dParams.domainmin.x) 
      + numprintf(4,0,rtime2dParams.domainmax.x)
      + numprintf(4,3,rtime2dParams.intergrid.x)
      + "_lat"       + numprintf(4,0,rtime2dParams.domainmin.y)
      + numprintf(4,0,rtime2dParams.domainmax.y)
      + numprintf(4,3,rtime2dParams.intergrid.y)
      + "_ts"        + numprintf(3,2,rtime2dParams.intstep) 
      + "_t"         + numprintf(4,0,rtime2dParams.tau) 
      + "_d"         + numprintf(2,0,rtime2dParams.seeddate.tm_mday)
      + numprintf(2,0,rtime2dParams.seeddate.tm_mon+1)
      + numprintf(2,0,rtime2dParams.seeddate.tm_year);
  }

  string nfilertime2d = rawname + ".data"; 
  ofstream ofilertime2d(nfilertime2d.c_str());
#ifdef DEBUG
  cout << "Save rtime values in file: " << nfilertime2d <<endl;
#endif
  for(q=0; q<rtime.size(); q++) {
    ofilertime2d<<rtime[q]<<endl;
  }  
  ofilertime2d.close();

#ifdef DEBUG
  cout << "[Complete]" << endl;
#endif

  /********************************************************************************
   * WRITING RESULT IN VTK FILE
   ********************************************************************************/

  string vtkfilertime2d = rawname + ".vtk";

#ifdef DEBUG
  cout << "Save ftle field in vtk file: " << vtkfilertime2d <<endl;
#endif

  ofstream vtkfile(vtkfilertime2d.c_str());
  vtkfile<<"# vtk DataFile Version 3.0"<<endl;
  vtkfile<<"Finite size Lyapunov exponent 2D"<<endl; 
  vtkfile<<"ASCII"<<endl;
  vtkfile<<"DATASET STRUCTURED_POINTS"<< endl;
  vtkfile<<"DIMENSIONS "<< dimension.i <<" "<< dimension.j <<" "<< 1 <<endl;
  vtkfile<<"ORIGIN "<<grid[0].x<<" "<<grid[0].y<<" "<< 0.0 <<endl;
  vtkfile<<"SPACING "<<rtime2dParams.intergrid.x<<" "<<rtime2dParams.intergrid.y<<" "<< 0.0 <<endl;
  vtkfile<<"POINT_DATA "<<dimension.i*dimension.j<<endl;
  vtkfile<<"SCALARS rtime2d double"<<endl;
  vtkfile<<"LOOKUP_TABLE default"<<endl;
  for(q=0; q<rtime.size(); q++) {
      vtkfile<<rtime[q]<<endl;
  }
  vtkfile<<"SCALARS qflag int"<<endl;
  vtkfile<<"LOOKUP_TABLE default"<<endl;
  for(q=0; q<numgridpoints; q++) {
    vtkfile<<qflag[q]<<endl;
  }
  vtkfile.close();

#ifdef DEBUG
  cout << "[Complete]" <<endl;
#endif 
  return 0;
}

string numprintf(int ndigits, int ndecimals, double number) {
  ostringstream stream;// it needs to include <sstream>
  double factor=pow(10,ndecimals); // it needs to include <cmath>
  stream << fixed; //it needs to include <iostream>
  stream << setfill('0') << setw(ndigits);
  stream << setprecision(0) << internal<< (number*factor);
  return stream.str();
}
