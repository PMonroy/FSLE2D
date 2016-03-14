#include "vectorIJ.hpp"

vectorIJ::vectorIJ(void) {//Void constructor
	i=0;
	j=0;
}
vectorIJ::vectorIJ(int ii,int jj) {//Constructor
	i=ii;
	j=jj;
}
vectorIJ::vectorIJ(const vectorIJ &a) {//Copy vector constructor
	i=a.i;
	j=a.j;
}
void vectorIJ::operator=(vectorIJ a) {//Copy vector
	i=a.i;
	j=a.j;
}

vectorIJ &vectorIJ::operator+=(vectorIJ a) {//Addition increment
    i += a.i;
    j += a.j;
    return *this;
}
vectorIJ &vectorIJ::operator-=(vectorIJ a) {//Subtraction increment
    i -= a.i;
    j -= a.j;
    return *this;
}
vectorIJ &vectorIJ::operator*=(vectorIJ a) {//Product increment
    i *= a.i;
    j *= a.j;
    return *this;
}
vectorIJ &vectorIJ::operator/=(vectorIJ a) {//Division increment
    i /= a.i;
    j /= a.j;
    return *this;
}
vectorIJ &vectorIJ::operator*=(int scalar) {//Scalar product increment  
    i *= scalar;
    j *= scalar;
    return *this;
}

vectorIJ  vectorIJ::operator+(vectorIJ a) {//Addition
  return vectorIJ(*this)+=a;
}
vectorIJ  vectorIJ::operator-(vectorIJ a) {//Subtraction
  return vectorIJ(*this)-=a;
}
vectorIJ  vectorIJ::operator*(vectorIJ a) {//Product
  return vectorIJ(*this)*=a;
}
vectorIJ  vectorIJ::operator/(vectorIJ a) {//Division
  return vectorIJ(*this)/=a;
}
vectorIJ operator*( int scalar,vectorIJ a) {//Scalar left product
  return vectorIJ(a) *= scalar;
} 
vectorIJ operator*(vectorIJ a, int scalar) {//Scalar right product
  return vectorIJ(a) *= scalar;
} 

vectorIJ::~vectorIJ() {//Destructor 

}    

ostream &operator<<(ostream &out, vectorIJ a) { //Format output
        out<<a.i<<" "<<a.j;
        return out;
}
istream &operator>>(istream &in, vectorIJ &a) {//Format input
  in>>a.i;
  in>>a.j;
  return in;
}

bool operator==(vectorIJ a,vectorIJ b) {//Boolean operator ==
  return (a.i==b.i && a.j==b.j);
}
