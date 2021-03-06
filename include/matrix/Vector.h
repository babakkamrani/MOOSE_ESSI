/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.8 $
// $Date: 2003/04/02 22:02:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/matrix/Vector.h,v $

                                                                        
// File: ~/matrix/Vector.h
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the class definition for Vector.
// Vector is a concrete class implementing the vector abstraction.
//
// What: "@(#) Vector.h, revA"

#ifndef Vector_h
#define Vector_h 

// Boris Jeremic 17Nov2008
#include <iostream>
//using namespace std;


class ID;

#define VECTOR_VERY_LARGE_VALUE 1.0e200

class Matrix; 
class Message;
class SystemOfEqn;

#include <BJtensor.h> // cannot use class as Tensor is itself defined in Tensor.h


class Vector
{
  public:
    // constructors and destructor
    Vector();
    Vector(int);
    Vector(const Vector &);    
    Vector(double *data, int size);
    ~Vector();

    // utility methods
    int setData(double *newData, int size);
    int Assemble(const Vector &V, const ID &l, double fact = 1.0);
    double Norm(void) const;
    inline int Size(void) const;
    int resize(int newSize);
    inline void Zero(void);
    int Normalize(void);
    
    int addVector(double factThis, const Vector &other, double factOther);
    int addMatrixVector(double factThis, const Matrix &m, const Vector &v, double factOther); 
    int addMatrixTransposeVector(double factThis, const Matrix &m, const Vector &v, double factOther);

    
    // overloaded operators
    inline double operator()(int x) const;
    inline double &operator()(int x);
    double operator[](int x) const;  // these two operator do bounds checks
    double &operator[](int x);
    Vector operator()(const ID &rows) const;
    Vector &operator=(const Vector  &V);
    Vector &operator=(const BJtensor  &T);
    
    Vector &operator+=(double fact);
    Vector &operator-=(double fact);
    Vector &operator*=(double fact);
    Vector &operator/=(double fact); 

    Vector operator+(double fact) const;
    Vector operator-(double fact) const;
    Vector operator*(double fact) const;
    Vector operator/(double fact) const;
    
    Vector &operator+=(const Vector &V);
    Vector &operator-=(const Vector &V);
    
    Vector operator+(const Vector &V) const;
    Vector operator-(const Vector &V) const;
    double operator^(const Vector &V) const;
    Vector operator/(const Matrix &M) const;    

    // methods added by Remo
    int  Assemble(const Vector &V, int init_row, double fact = 1.0);
    int  Extract (const Vector &V, int init_row, double fact = 1.0); 
  
//     friend OPS_Stream &operator<<(OPS_Stream &s, const Vector &V);
    // friend istream &operator>>(istream &s, Vector &V);    
    friend Vector operator*(double a, const Vector &V);
    
    friend class Message;
    friend class SystemOfEqn;
    friend class Matrix;
    friend class TCP_SocketNoDelay;    
    friend class TCP_Socket;
    friend class UDP_Socket;
    friend class MPI_Channel;
    friend class MySqlDatastore;
    friend class MySqlDataRecorder;
    friend class BerkeleyDbDatastore;
    
  private:
    static double VECTOR_NOT_VALID_ENTRY;
    int sz;
    double *theData;
    int fromFree;
};


/********* INLINED VECTOR FUNCTIONS ***********/
inline int 
Vector::Size(void) const 
{
  return sz;
}


inline void
Vector::Zero(void){
  for (int i=0; i<sz; i++) theData[i] = 0.0;
}


inline double 
Vector::operator()(int x) const
{
//#ifdef _G3DEBUG
  // check if it is inside range [0,sz-1]
  if (x < 0 || x >= sz) {
      std::cerr << "Vector::(loc) - loc " << x << " outside range [0, " << sz-1 ;
      return VECTOR_NOT_VALID_ENTRY;
  }
//#endif

      return theData[x];
}


inline double &
Vector::operator()(int x)
{
//#ifdef _G3DEBUG
    // check if it is inside range [0,sz-1]
  if (x < 0 || x >= sz) {
      std::cerr << "Vector::(loc) - loc " << x << " outside range [0, " << sz-1 ;
      return VECTOR_NOT_VALID_ENTRY;
  }
//#endif
  
  return theData[x];
}


#endif

