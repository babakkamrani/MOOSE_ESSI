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
**                                                                    **
** Additions and changes by:                                          **
**   Boris Jeremic (@ucdavis.edu)                                     **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.17 $
// $Date: 2004/07/20 22:39:02 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/NDMaterial.cpp,v $

// File: ~/material/NDMaterial.C
//
// Written: MHS
// Created: Feb 2000
// Revision: A
//
// Description: This file contains the class implementation for NDMaterial.
//
// What: "@(#) NDMaterial.C, revA"

#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <stresst.h>
#include <straint.h>
#include "BJtensor.h"


Matrix NDMaterial::errMatrix(1,1);
Vector NDMaterial::errVector(1);
BJtensor NDMaterial::errTensor(2, def_dim_2, 0.0 );
stresstensor NDMaterial::errstresstensor;
straintensor NDMaterial::errstraintensor;

NDMaterial::NDMaterial(int tag, int classTag)
{

}

NDMaterial::NDMaterial()
{

}

NDMaterial::~NDMaterial()
{

}


double
NDMaterial::getRho(void)
{
  return 0.0;
}

// BJ added 19June2002
double NDMaterial::getE(void)
  {
    return 0.0;
  }

// BJ added 19June2002
double NDMaterial::getnu(void)
  {
    return 0.0;
  }

// BJ added 19June2002
double NDMaterial::getpsi(void)
  {
    return 0.0;
  }

const Vector &
NDMaterial::getCommittedStress(void)
{
  return this->getStress();
}

const Vector &
NDMaterial::getCommittedStrain(void)
{
  return this->getStrain();
}

// methods to set and retrieve state.
int
NDMaterial::setTrialStrain(const Vector &v)
{
   cerr << "NDMaterial::setTrialStrain -- subclass responsibility\n";
   return -1;
}

int
NDMaterial::setTrialStrain(const Vector &v, const Vector &r)
{
   cerr << "NDMaterial::setTrialStrain -- subclass responsibility\n";
   return -1;
}

int
NDMaterial::setTrialStrainIncr(const Vector &v)
{
   cerr << "NDMaterial::setTrialStrainIncr -- subclass responsibility\n";
   return -1;
}

int
NDMaterial::setTrialStrainIncr(const Vector &v, const Vector &r)
{
   cerr << "NDMaterial::setTrialStrainIncr -- subclass responsibility\n";
   return -1;
}

const Matrix &
NDMaterial::getTangent(void)
{
   cerr << "NDMaterial::getTangent -- subclass responsibility\n";
   return errMatrix;
}

const Vector &
NDMaterial::getStress(void)
{
   cerr << "NDMaterial::getStress -- subclass responsibility\n";
   return errVector;
}

const Vector &
NDMaterial::getStrain(void)
{
   cerr << "NDMaterial::getStrain -- subclass responsibility\n";
   return errVector;
}

int
NDMaterial::setTrialStrain(const BJtensor &v)
{
   cerr << "NDMaterial::setTrialStrainIncr -- subclass responsibility\n";
   return -1;
}

int
NDMaterial::setTrialStrain(const BJtensor &v, const BJtensor &r)
{
   cerr << "NDMaterial::setTrialStrainIncr -- subclass responsibility\n";
   return -1;
}

int
NDMaterial::setTrialStrainIncr(const BJtensor &v)
{
   cerr << "NDMaterial::setTrialStrainIncr -- subclass responsibility\n";
   return -1;
}

int
NDMaterial::setTrialStrainIncr(const BJtensor &v, const BJtensor &r)
{
   cerr << "NDMaterial::setTrialStrainIncr -- subclass responsibility\n";
   return -1;
}

//Zhao (zcheng@ucdavis.edu)
// added Sept 22 2003 for Large Deformation, F is the Deformation Grandient
int
NDMaterial::setTrialF(const straintensor &f)
{
   cerr << "NDMaterial::setTrialF -- subclass responsibility\n";
   return -1;
}

int
NDMaterial::setTrialFIncr(const straintensor &df)
{
   cerr << "NDMaterial::setTrialF -- subclass responsibility\n";
   return -1;
}

int
NDMaterial::setTrialC(const straintensor &c)
{
   cerr << "NDMaterial::setTrialC -- subclass responsibility\n";
   return -1;
}

int
NDMaterial::setTrialCIncr(const straintensor &c)
{
   cerr << "NDMaterial::setTrialC -- subclass responsibility\n";
   return -1;
}

const stresstensor& NDMaterial::getPK1StressTensor(void)
{
   cerr << "NDMaterial::getPK1StressBJtensor -- subclass responsibility\n";
   return errstresstensor;
}

const stresstensor& NDMaterial::getCauchyStressTensor(void)
{
   cerr << "NDMaterial::getCauchyStressBJtensor -- subclass responsibility\n";
   return errstresstensor;
}

const straintensor& NDMaterial::getF(void)
{
   cerr << "NDMaterial::getF -- subclass responsibility\n";
   return errstraintensor;
}

const straintensor& NDMaterial::getC(void)
{
   cerr << "NDMaterial::getF -- subclass responsibility\n";
   return errstraintensor;
}

const straintensor& NDMaterial::getFp(void)
{
   cerr << "NDMaterial::getFp -- subclass responsibility\n";
   return errstraintensor;
}
// Only For Large Deformation, END////////////////////////////

const BJtensor &
NDMaterial::getTangentTensor(void)
{
   cerr << "NDMaterial::getTangentBJtensor -- subclass responsibility\n";
   return errTensor;
}

const stresstensor &NDMaterial::getStressTensor(void)
{
   cerr << "NDMaterial::getStressBJtensor -- subclass responsibility\n";
   return errstresstensor;
}

const straintensor &NDMaterial::getStrainTensor(void)
{
   cerr << "NDMaterial::getStrainBJtensor -- subclass responsibility\n";
   return errstraintensor;
}

const straintensor &NDMaterial::getPlasticStrainTensor(void)
{
   cerr << "NDMaterial::getPlasticStrainBJtensor -- subclass responsibility\n";
   return errstraintensor;
}


// Nima Tafazzoli added for checking the mesh (September 2012)
// int
// NDMaterial::CheckMesh(ofstream& checkmesh_file)
// {
//     return 0;
// }


