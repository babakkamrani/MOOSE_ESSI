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
// DESIGNER:          Mahdi Taiebat, Boris Jeremic
// PROGRAMMER:        Mahdi Taiebat, Boris Jeremic
// Note:
// DATE:              Spring 2007
// UPDATE HISTORY:
//
///////////////////////////////////////////////////////////////////////////////
//

// This is based on:
// G = G0*Pat*[(2.97-e)^2/(1+e)]*(p/Pat)^0.5
// K = K0*Pat*[(1.0+e)/e]*(p/Pat)^(0.66666666667)
//
// Constants:
// 1: G0:    Reference elastic shear modulus (same unit as stress);
// 2: K0:    Reference elastic bulk modulus (same unit as stress);
// 3: Pat:   Atmospheric pressure;
// 4: k_c;   cut-off factor,
//           for p < k_c*Pat, use p = k_c*Pat for calculation of G;
//           (a default value of k_c = 0.01 should work fine)
// 5: e0;    initial void ratio;


// Note-1: possible improvement: by default we can use K_c=0.01 i.e.
//         for p < 0.01*Pat we'll use p = 0.01*Pat for calculation
//         of the elastic bulk modulus.

#ifndef SANISAND_Elastic_CPP
#define SANISAND_Elastic_CPP

#include "SANISAND_Elastic.h"
#include <iostream>
using namespace std;


BJtensor SANISAND_Elastic::ElasticStiffness(4, def_dim_4, 0.0);

SANISAND_Elastic::SANISAND_Elastic(int G0_in,
                                   int K0_in,
                                   int Pat_in,
                                   int k_c_in,
                                   int e0_in,
                                   const stresstensor& initialStress,
                                   const straintensor& initialStrain)
    : ElasticState(ES_TAG_SANISAND, initialStress, initialStrain),
      G0_index(G0_in),
      K0_index(K0_in),
      Pat_index(Pat_in),
      k_c_index(k_c_in),
      e0_index(e0_in)
{

}

// Create a new
ElasticState* SANISAND_Elastic::newObj()
{
    ElasticState* Els = new SANISAND_Elastic(this->G0_index,
            this->K0_index,
            this->Pat_index,
            this->k_c_index,
            this->e0_index,
            this->Stress,
            this->Strain);
    return Els;
}
////////////////////////////////////////////////////////////////
SANISAND_Elastic::~SANISAND_Elastic()
{

}
////////////////////////////////////////////////////////////////

// Get Stiffness Tensor
const BJtensor& SANISAND_Elastic::getElasticStiffness(const MaterialParameter& MaterialParameter_in) const
{
    // Kronecker delta tensor
    BJtensor I2("I", 2, def_dim_2);

    BJtensor I_ijkl = I2("ij") * I2("kl");
    I_ijkl.null_indices();
    BJtensor I_ikjl = I_ijkl.transpose_1234_1324();
    BJtensor I_iljk = I_ijkl.transpose_1234_1423();
    BJtensor I4s = (I_ikjl + I_iljk) * 0.5;

    double G0 = getG0(MaterialParameter_in);
    double K0 = getK0(MaterialParameter_in);
    double Pat = getPat(MaterialParameter_in);
    double k_c = getk_c(MaterialParameter_in);
    double e0 = gete0(MaterialParameter_in);

    double epsilon_v = this->getStrain().Iinvariant1();
    double e = e0 + (1.0 + e0) * epsilon_v;
    double p_cal = this->getStress().p_hydrostatic();

    if (p_cal < 0.0)
    {
        p_cal = 0.0;
    }

    double p_cut = Pat * k_c;
    double p = (p_cal > p_cut) ? p_cal : p_cut;
    double G = G0 * Pat * ((2.97 - e) * (2.97 - e) / (1.0 + e)) * sqrt(p / Pat);
    double K = K0 * Pat * ((1.0 + e) / e) * pow(p / Pat, 0.6666666667);

    // Building elasticity tensor
    SANISAND_Elastic::ElasticStiffness = I_ijkl * (K - 2.0 * G / 3.0) + I4s * (2.0 * G);

    return SANISAND_Elastic::ElasticStiffness;
}

////////////////////////////////////////////////////////////////
const stresstensor& SANISAND_Elastic::getStress() const
{
    return ElasticState::Stress;
}

// Get G0
double SANISAND_Elastic::getG0(const MaterialParameter& MaterialParameter_in) const
{
    return MaterialParameter_in.getMaterial_Constant(G0_index - 1);
}

// Get K0
double SANISAND_Elastic::getK0(const MaterialParameter& MaterialParameter_in) const
{
    return MaterialParameter_in.getMaterial_Constant(K0_index - 1);
}

// Get Pat
double SANISAND_Elastic::getPat(const MaterialParameter& MaterialParameter_in) const
{
    return MaterialParameter_in.getMaterial_Constant(Pat_index - 1);
}

// Get k_cut
double SANISAND_Elastic::getk_c(const MaterialParameter& MaterialParameter_in) const
{
    return MaterialParameter_in.getMaterial_Constant(k_c_index - 1);
}

// Get e0
double SANISAND_Elastic::gete0(const MaterialParameter& MaterialParameter_in) const
{
    return MaterialParameter_in.getMaterial_Constant(e0_index - 1);
}


#endif

