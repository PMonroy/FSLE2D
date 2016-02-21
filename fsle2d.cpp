#include <cstdio>
#include <cstdlib>
#include <iostream> // cout, endl...  
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>
#include <iomanip> // setfill, setw

using namespace std;

#include "readparameters.h"
#include "gridconstruction.h" 
#include "constants.h"
#include "vflow.h" 
#include "integration.h"

// Maybe one day I use this:
//int (*velocity)(double ,vectorXY , vectorXY* ); 

string numprintf(int ndigits, int ndecimals, double number);

int main(int argc, char **argv)
{

  /**********************************************
   * READ COMAND LINE PARAMETERS
   **********************************************/
  string fnameparams;
  if(GetcmdlineParameters(argc, argv, &fnameparams))
    return 1;

  /* VERBOSE */
  if(verbose)
    {
      cout << "READ COMMAND LINE ARGUMENTS:" <<endl;
      cout << " Parameters file: " << fnameparams <<endl;
      cout << endl;
    }


  /**********************************************
   * READ PARAMETERS FROM FILE
   **********************************************/
  if(GetfileParameters(fnameparams))
    return 1;

  /* VERBOSE */
  if(verbose == 1)
    {
      cout << "PARAMETERS FROM FILE: "<< fnameparams <<endl; 
      cout << " vfield = "<<vfield<<endl ;
      cout << " domainmin = "<< domainmin<<endl;
      cout << " intergrid = " <<intergrid<<endl;
      cout << " domainmax = "<< domainmax<<endl;
      cout << " seeddate = "<< seeddate.tm_mday<<"-"<<seeddate.tm_mon+1<<"-"<<seeddate.tm_year<<endl ;
      cout << " intstep = "<<intstep<<endl ;
      cout << " tau = "<<tau<<endl ;
      cout << " deltamax = "<<deltamax<<endl ;
      cout << endl;
    }

  /**********************************************
   * GRID CONSTRUCTION
   **********************************************/
  vector<vectorXY> grid;
  int ni,nj;

  if(gridfsle2d(&grid, &ni, &nj, domainmin, intergrid, domainmax))
    {
      cout << "*****Error in contruction grid" << endl;
      return 1;
    }

  if(verbose == 1)
    {
      cout << "GRID CONSTRUCTION: "<<endl; 
      cout << " num. nodes = "<< grid.size()<<endl;
      cout << " ni = "<<ni<<endl;
      cout << " nj = "<<nj<<endl;
      cout << endl;
    }

  vector<int> neighbor;
  neighbor.reserve(4*grid.size());
  neighbor = neighbors(ni, nj);

  /**********************************************
   * COMPUTE INITIAL SEPARATION
   **********************************************/

  vectorXY delta,scalefactor;
  vector<double> ilength;
  unsigned int p;

  ilength.reserve(4*grid.size());

  for(unsigned int q=0; q<grid.size(); q++)
    {
      p=4*q;
      for(int dir=0; dir<4; dir++)
	{
	  if(neighbor[p+dir]!=-1)
	    {
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

  /**********************************************
   * INITIALIZE VARIBLES
   **********************************************/

  vector<double> exit_time;
  vector<double> response;
  vector<int> qcore;
  vector<int> qflag;

  //Reserve memory to improve efficiency
  exit_time.reserve(grid.size());
  response.reserve(grid.size());
  qcore.reserve(grid.size());
  qflag.reserve(grid.size());

  for(unsigned int q=0; q<grid.size(); q++)
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

  /**********************************************
   * SETUP TIME PARAMETERS
   **********************************************/

  struct tm *inidate = {0};

  double tstart;
  double tend;
  double h;

  int ntime = abs(tau);
  int ascnd = tau > 0;
  time_t seedtime = mktime(&seeddate); // Convert date to time in seconds (UTC) 

  if(ascnd)
    {
      time_t initime = seedtime;
      inidate = gmtime(&initime);
      tstart = 0.0;
      tend = (double) ntime;
      h=intstep;
    }
  else
    {
      time_t initime = seedtime - ntime*secondsday;
      inidate = gmtime(&initime); 
      tend = 0.0;
      tstart = (double) ntime;
      h=-1.0*intstep;
    }


  /**********************************************
   * LOAD VELOCITY FIELD
   **********************************************/
  if(verbose == 1) cout << "Loading grid:" << endl;

  if(loadvgrid(seeddate,vfield)!=0)  // Load velocity grid
    {
      cout << "***** Error reading grid" << endl;
      return 1;
    }

  if(verbose == 1) cout << "[Complete]" << endl;

  if(verbose == 1) cout << "Loading velocities:" << endl;

  if(loadvflow(*inidate, ntime+2, vfield)!=0)  // Load velocity field
    {
      cout << "***** Error reading velocities"<< endl;
      return 1;
    }
  if(verbose == 1) cout << "[Complete]"<<endl;

  /**********************************************
   * TIME LOOP
   **********************************************/
  
  vector<vectorXY> tracer = grid;  
  vector<double> length = ilength;

  double lengthmax;
  int dirmax;
  
  int q0, q1, q2, q3; 
  int qdir;

  unsigned int count;
  double t;

  if(verbose == 1) cout << "Calculation of response and exit_time:" << endl;

  for(t=tstart, count=0; ascnd==1?(t<tend):(t>=tend); t+=h,count++)
     {

       if(verbose == 1) cout << "count=" << count << "(" <<(tend-tstart)/h <<")"<<endl;

       /* Compute the position of tracer in time t=t+h*/
       for (unsigned int q = 0; q<tracer.size() ; q++)
	 {
	   //index of four q-neighboirs
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
	      (q3!=-1 && qflag[q3]==1))
	     {
	       if(RK4(t, h, &tracer[q], getvflow, vfield)==1)
		 {
		   qflag[q]=-1; // q is a non-integrable grid point
		   if(q0!=-1 && qflag[q0]!=-1) qflag[q0]=0;// We keep non-integrable neighbour grid points
		   if(q1!=-1 && qflag[q1]!=-1) qflag[q1]=0;
		   if(q2!=-1 && qflag[q2]!=-1) qflag[q2]=0;
		   if(q3!=-1 && qflag[q3]!=-1) qflag[q3]=0;
		 }
	     }
	 }

       /* Compute the relative distances */
       for(unsigned int q=0; q<tracer.size(); q++)
	 {
	   p=4*q;
	   for(int dir=0; dir<2; dir++)
	     {
	       qdir = neighbor[p+dir];
	       if(qdir==-1) // if grid point hasn't neighbour in direction dir,
		 continue;  // jump to the next step in the loop

	       if(qflag[q]==1 || qflag[qdir]==1)
		 {
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
       for(unsigned int q=0; q<tracer.size(); q++)
	 {
	   if(qflag[q]==1)
	     {
	       p=4*q;
	       lengthmax=length[p];// May we initialize lengthmax before?
	       for(int dir=1; dir<4; dir++)
		 {
		   if(length[p+dir]>lengthmax)
		     {
		       lengthmax=length[p+dir];
		       dirmax=dir;
		     }
		 }

	       if(lengthmax>deltamax)
		 {
		   qflag[q]=0;
		   exit_time[q]=(count+1.0)*h;
		   response[q]=lengthmax/ilength[p+dirmax];
		 }
	     }
	 }
     }
  if(verbose==1)  cout << "[Complete]" << endl;

  /**********************************************************
   * COMPUTE FINITE SIZE LYAPUNOV EXPONENT
   **********************************************************/

  vector<double> fsle;
  fsle.reserve(grid.size());

  for(unsigned int q=0; q<qcore.size(); q++)
    {
      fsle.push_back(log(response[qcore[q]])/exit_time[qcore[q]]);
    }

  // Save fsle 2d grid in a file
  string nfilegridfsle2d = 
    "fsle2d-dmin" + numprintf(4,0,domainmin.x) +
    "_"            + numprintf(4,0,domainmin.y) +
    "_dmax"       + numprintf(4,0,domainmax.x) +
    "_"            + numprintf(4,0,domainmax.y) +
    "_res"        + numprintf(4,3,intergrid.x) +
    "_"            + numprintf(4,3,intergrid.y) +
    ".grid";
  ofstream ofilegridfsle2d(nfilegridfsle2d.c_str());

  if(verbose==1) cout << "Save grid in file: " << nfilegridfsle2d <<endl;

  for(unsigned int q=0; q<qcore.size(); q++) 
    {
      ofilegridfsle2d<<grid[q].x<<" "<<grid[q].y<<endl;
    }
  ofilegridfsle2d.close();

  if(verbose==1) cout << "[Complete]" << endl;

  // save fsle 2d values in a file
  string nfilefsle2d = 
    "fsle2d_vm"   + numprintf(1,0,vfield) +
    "_date"       + numprintf(2,0,seeddate.tm_year) + 
    "-"           + numprintf(2,0,seeddate.tm_mon+1)+
    "-"           + numprintf(2,0,seeddate.tm_mday) +
    "_rlon"       + numprintf(4,0,domainmin.x) +
    "_"           + numprintf(4,0,domainmax.x) +
    "_rlat"       + numprintf(4,0,domainmin.y) +
    "_"           + numprintf(4,0,domainmax.x) +
    "_res"        + numprintf(4,3,intergrid.x) +
    "_"           + numprintf(4,3,intergrid.y) +
    "_tau"        + numprintf(4,0,tau) +
    "_h"          + numprintf(4,3,intstep) +
    "_dmax"       + numprintf(3,0, deltamax/1000.0) +
    ".data";

  if(verbose==1) cout << "Save fsle values in file: " << nfilefsle2d <<endl;

  ofstream ofilefsle2d(nfilefsle2d.c_str());
  for(unsigned int q=0; q<fsle.size(); q++) 
    {
      ofilefsle2d<<fsle[q]<<endl;
    }  
  ofilefsle2d.close();

  if(verbose==1)  cout << "[Complete]" << endl;


 /**********************************************************
   * WRITING RESULT IN VTK FILE
   **********************************************************/

  string vtkfilefsle2d = 
    "fsle2d_vm"   + numprintf(1,0,vfield) +
    "_date"       + numprintf(2,0,seeddate.tm_year) + 
    "-"           + numprintf(2,0,seeddate.tm_mon+1)+
    "-"           + numprintf(2,0,seeddate.tm_mday) +
    "_rlon"       + numprintf(4,0,domainmin.x) +
    "_"           + numprintf(4,0,domainmax.x) +
    "_rlat"       + numprintf(4,0,domainmin.y) +
    "_"           + numprintf(4,0,domainmax.x) +
    "_res"        + numprintf(4,3,intergrid.x) +
    "_"           + numprintf(4,3,intergrid.y) +
    "_tau"        + numprintf(4,0,tau) +
    "_h"          + numprintf(4,3,intstep) +
    "_dmax"       + numprintf(3,0, deltamax/1000.0) +
    ".vtk";

  if(verbose==1)  cout << "Save ftle field in vtk file: " << vtkfilefsle2d <<endl;

  /*ofstream vtkfile(vtkfilefsle2d.c_str());
  vtkfile<<"# vtk DataFile Version 3.0"<<endl;
  vtkfile<<"Finite size Lyapunov exponent 2D"<<endl; 
  vtkfile<<"ASCII"<<endl;
  vtkfile<<"DATASET STRUCTURED_POINTS"<< endl;
  vtkfile<<"DIMENSIONS "<< ni-2 <<" "<< nj-2 <<" "<< 1 <<endl;
  vtkfile<<"ORIGIN "<<grid[qcore[0]].x<<" "<<grid[qcore[0]].y<<" "<< 0.0 <<endl;
  vtkfile<<"SPACING "<<intergrid.x<<" "<<intergrid.y<<" "<< 0.0 <<endl;
  vtkfile<<"POINT_DATA "<<(ni-2)*(nj-2)<<endl;
  vtkfile<<"SCALARS fsle2d double"<<endl;
  vtkfile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<fsle.size(); q++) 
    {
      vtkfile<<fsle[q]<<endl;
    }
  vtkfile<<"SCALARS qflag int"<<endl;
  vtkfile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<qcore.size(); q++) 
    {
      vtkfile<<qflag[qcore[q]]<<endl;
    }
    vtkfile.close();*/

  ofstream vtkfile(vtkfilefsle2d.c_str());
  vtkfile<<"# vtk DataFile Version 3.0"<<endl;
  vtkfile<<"Finite size Lyapunov exponent 2D"<<endl; 
  vtkfile<<"ASCII"<<endl;
  vtkfile<<"DATASET STRUCTURED_POINTS"<< endl;
  vtkfile<<"DIMENSIONS "<< ni-1 <<" "<< nj-1 <<" "<< 1 <<endl;
  vtkfile<<"ORIGIN "<<grid[qcore[0]].x-(intergrid.x/2.0)<<" "<<grid[qcore[0]].y-(intergrid.y/2.0)<<" "<< 0.0 <<endl;
  vtkfile<<"SPACING "<<intergrid.x<<" "<<intergrid.y<<" "<< 0.0 <<endl;
  vtkfile<<"CELL_DATA "<<(ni-2)*(nj-2)<<endl;
  vtkfile<<"SCALARS fsle2d double"<<endl;
  vtkfile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<fsle.size(); q++) 
    {
      vtkfile<<fsle[q]<<endl;
    }
  vtkfile<<"SCALARS qflag int"<<endl;
  vtkfile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<qcore.size(); q++) 
    {
      vtkfile<<qflag[qcore[q]]<<endl;
    }
    vtkfile.close();

  if(verbose==1)  cout << "[Complete]" <<endl;

  return 0;
}

string numprintf(int ndigits, int ndecimals, double number)
{
  ostringstream stream;// it needs to include <sstream>
  double factor=pow(10,ndecimals); // it needs to include <cmath>
  stream << fixed; //it needs to include <iostream>
  stream << setfill('0') << setw(ndigits);
  stream << setprecision(0) << internal<< (number*factor);
  return stream.str();
}
