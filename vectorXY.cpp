#include "vectorXY.hpp"

vectorXY::vectorXY(void) { // Void Constructor
  x=0.0;
  y=0.0;
}
vectorXY::vectorXY(double xi,double yi) { // Constructor
  x=xi;
  y=yi;
}
vectorXY::vectorXY(const vectorXY &a){ // Copy Vector Constructor
  x=a.x;
  y=a.y;
}
void vectorXY::operator=(vectorXY a){ // Copy Vector
  x=a.x;
  y=a.y;
}
vectorXY &vectorXY::operator+=(vectorXY a)
{
  x += a.x;
  y += a.y;
  return *this;
}

vectorXY &vectorXY::operator-=(vectorXY a)
{
  x -= a.x;
  y -= a.y;
  return *this;
}

vectorXY &vectorXY::operator*=(double scalar)
{
  x *= scalar;
  y *= scalar;
  return *this;
}

vectorXY &vectorXY::operator*=(vectorXY a)
{
  x *= a.x;
  y *= a.y;
  return *this;
}

vectorXY &vectorXY::operator/=(vectorXY a)
{
  x /= a.x;
  y /= a.y;
  return *this;
}

vectorXY vectorXY::operator+(vectorXY a) // Addition
{
  return vectorXY(*this)+=a;
}

vectorXY vectorXY::operator-(vectorXY a) // Subtraction
{
  return vectorXY(*this)-=a;
}

vectorXY vectorXY::operator*(vectorXY a) // Product
{
  return vectorXY(*this)*=a;
}

vectorXY vectorXY::operator/(vectorXY a) // Division
{
  return vectorXY(*this)/=a;
}

vectorXY::~vectorXY(){} // Destructor

vectorXY operator*( double scalar,vectorXY a) 
{
  return vectorXY(a) *= scalar;
}

vectorXY operator*(vectorXY a, double scalar) 
{
  return vectorXY(a) *= scalar;
}

ostream &operator<<(ostream &out, vectorXY a) //output
{
  out<<a.x<<" "<<a.y;
  return out;
}

istream &operator>>(istream &in, vectorXY &a) //input
{
  in>>a.x;
  in>>a.y;
  return in;
}
bool operator==(vectorXY a,vectorXY b)
{
  return (a.x==b.x && a.y==b.y);
}
