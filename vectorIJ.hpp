#ifndef VECTORIJ
#define VECTORIJ

#include <iostream>

using namespace std;

class  vectorIJ {// Class int 3D vector

public:
  int i,j; // Order list of 2 elements |i|j|

  vectorIJ(void);                 //Void Constructor
  vectorIJ(int ii,int jj); //Constructor
  vectorIJ(const vectorIJ &a);   //Copy vector constructor
  void operator=(vectorIJ a);     //Copy vector

  vectorIJ &operator+=(vectorIJ a); //Addition increment
  vectorIJ &operator-=(vectorIJ a); //Subtraction increment
  vectorIJ &operator*=(vectorIJ a); //Product increment
  vectorIJ &operator/=(vectorIJ a); //Division increment
  vectorIJ &operator*=(int scalar);  //Scalar product increment

  vectorIJ operator+( vectorIJ a);  // Addition
  vectorIJ operator-( vectorIJ a);  // Subtraction
  vectorIJ operator*( vectorIJ a);  // Product
  vectorIJ operator/( vectorIJ a);  // Division

  friend ostream &operator<<(ostream &out, vectorIJ a); //Format output
  friend istream &operator>>(istream &in, vectorIJ &a); //Format input

  ~vectorIJ();	                            // Destructor
};

vectorIJ  operator*(int scalar, vectorIJ a); // Scalar left product
vectorIJ  operator*(vectorIJ a,int scalar);  // Scalar right product

bool operator==(vectorIJ a,vectorIJ b); //Boolean operator ==

#endif 
