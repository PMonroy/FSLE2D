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
struct fsle2dParameters {

  const vectorXY domainmin;
  const vectorXY domainmax;
  const vectorXY intergrid;
  const struct tm seeddate;
  const double intstep;
  const int tau;
  const double deltamax;

  // Define a constructor that will load stuff from a configuration file.
  fsle2dParameters(const string & fsle2dParamsFileName)
  :domainmin(getVectorXYParam(fsle2dParamsFileName, "domainmin"))
  ,domainmax(getVectorXYParam(fsle2dParamsFileName, "domainmax"))
  ,intergrid(getVectorXYParam(fsle2dParamsFileName, "intergrid")) 
  ,seeddate(getDateParam(fsle2dParamsFileName, "seeddate")) 
  ,intstep(getDoubleParam(fsle2dParamsFileName, "intstep")) 
  ,tau(getIntParam(fsle2dParamsFileName, "tau"))
  ,deltamax(getDoubleParam(fsle2dParamsFileName, "deltamax"))
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

  const fsle2dParameters fsle2dParams(fnameparams);

  /* VERBOSE */
#ifdef DEBUG
    cout << "FSLE2D parameters from file: ";
    cout << fnameparams <<endl; 
    cout << " domainmin = "<< fsle2dParams.domainmin<<endl;
    cout << " intergrid = " <<fsle2dParams.intergrid<<endl;
    cout << " domainmax = "<< fsle2dParams.domainmax<<endl;
    cout << " seeddate = "<< fsle2dParams.seeddate.tm_mday;
    cout << " " << fsle2dParams.seeddate.tm_mon+1;
    cout << " " << fsle2dParams.seeddate.tm_year<<endl;
    cout << " intstep = "<<fsle2dParams.intstep<<endl ;
    cout << " tau = "<<fsle2dParams.tau<<endl ;
    cout << " deltamax = "<<fsle2dParams.deltamax<<endl ;
    cout << endl;
#endif

  /********************************************************************************
   * GRID CONSTRUCTION
   ********************************************************************************/

  vector<vectorXY> grid;
  vectorIJ dimension;
  vectorXY griddomainmin=fsle2dParams.domainmin;
  vectorXY griddomainmax=fsle2dParams.domainmax;

  griddomainmin-=fsle2dParams.intergrid; 
  griddomainmax+=fsle2dParams.intergrid;

  if(MakeRegularGrid(&grid, &dimension, griddomainmin, fsle2dParams.intergrid, griddomainmax)) {
    cout << "Error in contruction grid" << endl;
    return 1;
  }
  unsigned int numgridpoints=grid.size();

  vector<int> neighbor;
  neighbor.reserve(4*numgridpoints);
  neighbor = GetNeighbors(dimension);
  unsigned int numneighbors=neighbor.size();
 

#ifdef DEBUG
  cout << "GRID CONSTRUCTION: "<<endl; 
  cout << " num. nodes = "<< grid.size()<<endl;
  cout << " dim.i = "<<dimension.i<<endl;
  cout << " dim.j = "<<dimension.j<<endl;
  cout << endl;
#endif
  /********************************************************************************
   * COMPUTE INITIAL SEPARATION
   ********************************************************************************/

  vectorXY delta,scalefactor;
  vector<double> ilength;
  unsigned int p;
  unsigned int q;
  unsigned int dir;

  ilength.reserve(numneighbors);

  for(q=0; q<numgridpoints; q++) {
    p=4*q;
    for(dir=0; dir<4; dir++) {
      if(neighbor[p+dir]!=-1) {
	delta=grid[neighbor[p+dir]]-grid[q];
	
	delta.x=rads*delta.x;
	delta.y=rads*delta.y;
	
	scalefactor.x=rearth*cos(rads*grid[q].y); 
	scalefactor.y=rearth; 
	
	delta=delta*scalefactor;
	delta*=delta;
	
	ilength.push_back(sqrt(delta.x+delta.y));
      }
      else
	ilength.push_back(-1.0);
    }
  }

  /********************************************************************************
   * INITIALIZE VARIABLES
   ********************************************************************************/

  vector<double> exit_time;
  vector<double> response;
  vector<int> qcore;
  vector<int> qflag;

  //Reserve memory to improve efficiency
  exit_time.reserve(numgridpoints);
  response.reserve(numgridpoints);
  qcore.reserve(numgridpoints);
  qflag.reserve(numgridpoints);

  for(q=0; q<numgridpoints; q++)
    {
      p=4*q;
      if(neighbor[p]!=-1 && neighbor[p+1]!=-1 && neighbor[p+2]!=-1 && neighbor[p+3]!=-1)
	{
	  qflag.push_back(1);// Initially we can compute fsle in this points (4=num. neighbours)
	  qcore.push_back(q);// qcore contains all the index of the later points
	}
      else
	{
	  qflag.push_back(0); // Border points: we can not compute fsle on it (4>num. neighbours)
	}

      response.push_back(1.0);
      exit_time.push_back(1.0);
    }

  /********************************************************************************
   * SETUP TIME PARAMETERS
   ********************************************************************************/

  struct tm *inidate = {0};

  double tstart;
  double tend;
  double h;

  int ntime = abs(fsle2dParams.tau);
  int ascnd = fsle2dParams.tau > 0;
  struct tm datebuff = fsle2dParams.seeddate; // Buffer var to avoid error "cast const tm* to tm*"
  time_t seedtime = mktime(&datebuff); // Convert date to time in seconds (UTC) 

  if(ascnd){
    time_t initime = seedtime;
    inidate = gmtime(&initime);
    tstart = 0.0;
    tend = (double) ntime;
    h= fsle2dParams.intstep;
  }
  else{
    time_t initime = seedtime - ntime*secondsday;
    inidate = gmtime(&initime); 
    tend = 0.0;
    tstart = (double) ntime;
    h=-1.0*fsle2dParams.intstep;
  }

  /********************************************************************************
   * LOAD VELOCITY FIELD
   ********************************************************************************/

#ifdef DEBUG
  cout << "Loading velocity grid:" << endl;
#endif

 vectorXY meanvel(0.2,0.2);
 if(LoadVGrid(fsle2dParams.seeddate,
	      fsle2dParams.domainmin,
	      fsle2dParams.domainmax, 
	      meanvel, 
	      abs(fsle2dParams.tau))!=0){//Load velocity grid
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
  
  vector<vectorXY> tracer = grid;  
  vector<double> length = ilength;

  double lengthmax;
  int dirmax;
  
  int q0, q1, q2, q3; 
  int qdir;

  unsigned int count;
  double t;

#ifdef DEBUG
  cout << "Calculation of response and exit_time:" << endl;
#endif

  for(t=tstart, count=0; ascnd==1?(t<tend):(t>=tend); t+=h,count++) {

#ifdef DEBUG
    cout << "count=" << count << "(" <<(tend-tstart)/h <<")"<<endl;
#endif

       /* Compute the position of tracer in time t=t+h*/
    for (q = 0; q<numgridpoints ; q++) {
      //index of four q-neighbors
      q0 = neighbor[4*q];
      q1 = neighbor[4*q+1];
      q2 = neighbor[4*q+2];
      q3 = neighbor[4*q+3];
      
      if(qflag[q]==-1)
	continue;// skip the rest of the loop if q is a non-integrable point 
      
      if((qflag[q]==1) ||
	 (q0!=-1 && qflag[q0]==1) ||
	 (q1!=-1 && qflag[q1]==1) ||
	 (q2!=-1 && qflag[q2]==1) ||
	 (q3!=-1 && qflag[q3]==1)) {
	  if(RK4(t, h, &tracer[q], GetVFlow)==1) {
	    qflag[q]=-1; // q is a non-integrable grid point
	    if(q0!=-1 && qflag[q0]!=-1) qflag[q0]=0;// We keep non-integrable neighbour grid points
	    if(q1!=-1 && qflag[q1]!=-1) qflag[q1]=0;
	    if(q2!=-1 && qflag[q2]!=-1) qflag[q2]=0;
	    if(q3!=-1 && qflag[q3]!=-1) qflag[q3]=0;
	  }
      }
    }
    
    /* Compute the relative distances */
    for(q=0; q<numgridpoints; q++) {
      p=4*q;
      for(int dir=0; dir<2; dir++) {
	qdir = neighbor[p+dir];
	if(qdir==-1) // if grid point hasn't neighbour in direction dir,
	  continue;  // jump to the next step in the loop
	
	if(qflag[q]==1 || qflag[qdir]==1) {
	  delta=tracer[qdir]-tracer[q];
	  
	  delta.x=rads*delta.x;
	  delta.y=rads*delta.y;
	  
	  scalefactor.x=rearth*cos(rads*grid[q].y); 
	  scalefactor.y=rearth;
	  
	  delta=delta*scalefactor;
	  delta*=delta;
	  
	  length[p+dir]=sqrt(delta.x+delta.y);
	  length[4*qdir+2]=length[p+dir];
	}
      }
    }
    
    /* Compute the max length*/
    for(q=0; q<numgridpoints; q++) {
      if(qflag[q]==1) {
	p=4*q;
	lengthmax=length[p];// May we initialize lengthmax before?
	dirmax=0;
	for(int dir=1; dir<4; dir++) {
	  if(length[p+dir]>lengthmax) {
	    lengthmax=length[p+dir];
	    dirmax=dir;
	  }
	}
	
	if(lengthmax>fsle2dParams.deltamax)  {
	  qflag[q]=0;
	  exit_time[q]=(count+1.0)*h;
	  response[q]=lengthmax/ilength[p+dirmax];
	}
      }
    }
  }

  //free Vflow resources
  FreeVFlow((unsigned int)(ntime+2));

#ifdef DEBUG
  cout << "[Complete]" << endl;
#endif
 
  /********************************************************************************
   * COMPUTE FINITE SIZE LYAPUNOV EXPONENT
   ********************************************************************************/

  vector<double> fsle;
  fsle.reserve(numgridpoints);

  for(q=0; q<qcore.size(); q++){
    fsle.push_back(log(response[qcore[q]])/exit_time[qcore[q]]);
  }

  /********************************************************************************
   * SAVE RESULTS IN FILES
   ********************************************************************************/

  // Save fsle 2d grid in a file
  string nfilegridfsle2d =     
    "lon"    + numprintf(4,0,fsle2dParams.domainmin.x) 
    + numprintf(4,0,fsle2dParams.domainmax.x)
    + numprintf(4,3,fsle2dParams.intergrid.x)
    + "_lat"       + numprintf(4,0,fsle2dParams.domainmin.y)
    + numprintf(4,0,fsle2dParams.domainmax.y)
    + numprintf(4,3,fsle2dParams.intergrid.y)
    + ".grid";

#ifdef DEBUG
  cout << "Save grid in file: " << nfilegridfsle2d <<endl;
#endif

  if(!ifstream(nfilegridfsle2d.c_str())){ // Check if file exists
    ofstream ofilegridfsle2d(nfilegridfsle2d.c_str());
    for(q=0; q<qcore.size(); q++) {
      ofilegridfsle2d<<grid[qcore[q]].x<<" "<<grid[qcore[q]].y<<endl;
    }
    ofilegridfsle2d.close();
  }
  
#ifdef DEBUG
  cout << "[Complete]" << endl;
#endif

  // save fsle 2d values in a file
  string rawname;

  if(namefileflag==1){
    size_t lastdot = fnameparams.find_last_of(".");
    if(lastdot == string::npos){
      rawname="fsle2d_"+fnameparams;
    }
    else{
      rawname="fsle2d_"+fnameparams.substr(0,lastdot);
    }
  }
  else{
    rawname = "fsle2d_lon"    + numprintf(4,0,fsle2dParams.domainmin.x) 
      + numprintf(4,0,fsle2dParams.domainmax.x)
      + numprintf(4,3,fsle2dParams.intergrid.x)
      + "_lat"       + numprintf(4,0,fsle2dParams.domainmin.y)
      + numprintf(4,0,fsle2dParams.domainmax.y)
      + numprintf(4,3,fsle2dParams.intergrid.y)
      + "_dmax"        + numprintf(3,0,fsle2dParams.deltamax/1000.0)
      + "_ts"        + numprintf(3,2,fsle2dParams.intstep) 
      + "_t"         + numprintf(4,0,fsle2dParams.tau) 
      + "_d"         + numprintf(2,0,fsle2dParams.seeddate.tm_mday)
      + numprintf(2,0,fsle2dParams.seeddate.tm_mon+1)
      + numprintf(2,0,fsle2dParams.seeddate.tm_year);
  }

  string nfilefsle2d = rawname + ".data"; 
  ofstream ofilefsle2d(nfilefsle2d.c_str());
#ifdef DEBUG
  cout << "Save fsle values in file: " << nfilefsle2d <<endl;
#endif
  for(q=0; q<fsle.size(); q++) {
    ofilefsle2d<<fsle[q]<<endl;
  }  
  ofilefsle2d.close();

#ifdef DEBUG
  cout << "[Complete]" << endl;
#endif

  /********************************************************************************
   * WRITING RESULT IN VTK FILE
   ********************************************************************************/

  string vtkfilefsle2d = rawname + ".vtk";

#ifdef DEBUG
  cout << "Save ftle field in vtk file: " << vtkfilefsle2d <<endl;
#endif

  ofstream vtkfile(vtkfilefsle2d.c_str());
  vtkfile<<"# vtk DataFile Version 3.0"<<endl;
  vtkfile<<"Finite size Lyapunov exponent 2D"<<endl; 
  vtkfile<<"ASCII"<<endl;
  vtkfile<<"DATASET STRUCTURED_POINTS"<< endl;
  vtkfile<<"DIMENSIONS "<< dimension.i-2 <<" "<< dimension.j-2 <<" "<< 1 <<endl;
  vtkfile<<"ORIGIN "<<grid[qcore[0]].x<<" "<<grid[qcore[0]].y<<" "<< 0.0 <<endl;
  vtkfile<<"SPACING "<<fsle2dParams.intergrid.x<<" "<<fsle2dParams.intergrid.y<<" "<< 0.0 <<endl;
  vtkfile<<"POINT_DATA "<<(dimension.i-2)*(dimension.j-2)<<endl;
  vtkfile<<"SCALARS fsle2d double"<<endl;
  vtkfile<<"LOOKUP_TABLE default"<<endl;
  for(q=0; q<fsle.size(); q++) {
      vtkfile<<fsle[q]<<endl;
  }
  vtkfile<<"SCALARS qflag int"<<endl;
  vtkfile<<"LOOKUP_TABLE default"<<endl;
  for(q=0; q<qcore.size(); q++) {
    vtkfile<<qflag[qcore[q]]<<endl;
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
