#ifndef VECTORXY
#define VECTORXY
#include <iostream>

using namespace std;

class vectorXY
{
 public:
  double x,y; // Order list of 3 elements |x|y|
  vectorXY(void); // Zero Vector Constructor
  vectorXY(double xi,double yi); // Constructor
  vectorXY(const vectorXY &a); // Copy Vector Constructor
  void operator = (vectorXY a);
  vectorXY &operator += (vectorXY a);
  vectorXY &operator -= (vectorXY a);
  vectorXY &operator *= (vectorXY a);
  vectorXY &operator *= (double scalar);
  vectorXY &operator /= (vectorXY a);
  vectorXY operator + (vectorXY a); // Addition
  vectorXY operator - (vectorXY a); // Subtraction
  vectorXY operator * (vectorXY a); // Product element by element
  vectorXY operator / (vectorXY a); // Product element by element
  friend ostream &operator << (ostream &out, vectorXY a);
  friend istream &operator >> (istream &in, vectorXY &a);
  ~vectorXY(); // Destructor
};

vectorXY operator*(double scalar, vectorXY a); // Product element by element
vectorXY operator*(vectorXY a,double scalar); // Product element by element
bool operator==(vectorXY a,vectorXY b);
#endif
