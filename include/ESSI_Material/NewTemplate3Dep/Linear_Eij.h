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

#ifndef Linear_Eij_H
#define Linear_Eij_H

#include "TensorEvolution.h"
#include "classTags.h"

class Linear_Eij : public TensorEvolution
{
    public:

        Linear_Eij() : TensorEvolution(TE_TAG_Linear_Eij) {}; //Guanzhou added for parallel processing

        Linear_Eij(int LinearFactor_index_in);

        TensorEvolution* newObj();
        ~Linear_Eij();

        const straintensor& Hij(const PlasticFlow& plastic_flow, const stresstensor& Stre,
                                const straintensor& Stra, const MaterialParameter& material_parameter);

        const tensor& DHij_Ds(const PlasticFlow& plastic_flow, const stresstensor& Stre,
                              const straintensor& Stra, const MaterialParameter& material_parameter);

        const tensor& DHij_Diso(const PlasticFlow& plastic_flow, const stresstensor& Stre,
                                const straintensor& Stra, const MaterialParameter& material_parameter);

        const tensor& DHij_Dkin(const PlasticFlow& plastic_flow, const stresstensor& Stre,
                                const straintensor& Stra, const MaterialParameter& material_parameter);

    private:

        double getLinearFactor(const MaterialParameter& material_parameter) const;

    private:

        int LinearFactor_index;
};


#endif

