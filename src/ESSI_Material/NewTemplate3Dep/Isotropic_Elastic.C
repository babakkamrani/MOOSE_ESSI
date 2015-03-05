///////////////////////////////////////////////////////////////////////////////
//   COPYLEFT (C): Woody's viral GPL-like license (by BJ):
//                 ``This    source  code is Copyrighted in
//                 U.S.,  for  an  indefinite  period,  and anybody
//                 caught  using it without our permission, will be
//                 mighty good friends of ourn, cause we don't give
//                 a  darn.  Hack it. Compile it. Debug it. Run it.
//                 Yodel  it.  Enjoy it. We wrote it, that's all we
//                 wanted to do.''
//
//
// COPYRIGHT (C):     :-))
// PROJECT:           Object Oriented Finite Element Program
// FILE:
// CLASS:
// MEMBER FUNCTIONS:
//
// MEMBER VARIABLES
//
// PURPOSE:
//
// RETURN:
// VERSION:
// LANGUAGE:          C++
// TARGET OS:
// DESIGNER:          Zhao Cheng, Boris Jeremic
// PROGRAMMER:        Zhao Cheng,
// DATE:              Fall 2005
// UPDATE HISTORY:    Guanzhou Jie updated for parallel Dec 2006
//
///////////////////////////////////////////////////////////////////////////////
//

#ifndef Isotropic_Elastic_CPP
#define Isotropic_Elastic_CPP

#define ES_TAG_IsotropicElastic 3205

#include "Isotropic_Elastic.h"
#include <iostream>
using namespace std;


BJtensor Isotropic_Elastic::ElasticStiffness(4, def_dim_4, 0.0);

Isotropic_Elastic::Isotropic_Elastic(int E_in,
                                     int v_in,
                                     const stresstensor& initialStress,
                                     const straintensor& initialStrain)
    : ElasticState(ES_TAG_IsotropicElastic, initialStress, initialStrain),
      E_index(E_in),
      v_index(v_in)
{

}


// Create a new
ElasticState* Isotropic_Elastic::newObj()
{
    ElasticState* Els = new  Isotropic_Elastic(this->E_index,
            this->v_index,
            this->Stress,
            this->Strain);
    return Els;
}
////////////////////////////////////////////////////////////////
Isotropic_Elastic::~Isotropic_Elastic()
{

}
////////////////////////////////////////////////////////////////
// Get Stiffness Tensor
const BJtensor& Isotropic_Elastic::getElasticStiffness(const MaterialParameter& MaterialParameter_in) const
{
    // Kronecker delta tensor
    BJtensor I2("I", 2, def_dim_2);

    BJtensor I_ijkl = I2("ij") * I2("kl");
    I_ijkl.null_indices();
    BJtensor I_ikjl = I_ijkl.transpose0110();
    BJtensor I_iljk = I_ijkl.transpose0111();
    BJtensor I4s = (I_ikjl + I_iljk) * 0.5;

    double E = getE(MaterialParameter_in);
    double v = getv(MaterialParameter_in);

    if (E < 0.0 || v < -1.0 || v >= 0.5)
    {
        cerr << "Isotropic_Elastic: Invalid Input. \n";
        exit (1);
    }

    // Building elasticity tensor
    Isotropic_Elastic::ElasticStiffness = I_ijkl * ( E * v / ( (1.0 + v) * (1.0 - 2.0 * v) ) ) + I4s * ( E / (1.0 + v) );

    return Isotropic_Elastic::ElasticStiffness;
}

// Get Young's modulus
double Isotropic_Elastic::getE(const MaterialParameter& MaterialParameter_in) const
{
    if ( E_index > MaterialParameter_in.getNum_Material_Constant() || E_index < 1)
    {
        cerr << "Isotropic_Elastic: Invalid Input. \n ";
        exit (1);
    }
    else
    {
        return MaterialParameter_in.getMaterial_Constant(E_index - 1);
    }
}


// Get Poisson's ratio
double Isotropic_Elastic::getv(const MaterialParameter& MaterialParameter_in) const
{
    if ( v_index > MaterialParameter_in.getNum_Material_Constant() || v_index  < 1)
    {
        cerr << "Isotropic_Elastic: Invalid Input. \n ";
        exit (1);
    }
    else
    {
        return MaterialParameter_in.getMaterial_Constant(v_index - 1);
    }
}

#endif

