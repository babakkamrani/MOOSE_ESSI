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
#ifndef PlasticFlow_H
#define PlasticFlow_H

#include <stresst.h>
#include <straint.h>
#include "MaterialParameter.h"
#include "classTags.h"


class PlasticFlow
{
    public:

        PlasticFlow(int clsTag);

        virtual ~PlasticFlow();

        virtual PlasticFlow* newObj() = 0;

        virtual const char* getPlasticFlowType(void) const = 0;

        virtual const straintensor& PlasticFlowTensor(const stresstensor& Stre,
                const straintensor& Stra,
                const MaterialParameter& MaterialParameter_in) const = 0;

        virtual const tensor& Dm_Ds(const stresstensor& Stre,
                                    const straintensor& Stra,
                                    const MaterialParameter& MaterialParameter_in) const;

        virtual const tensor& Dm_Diso(const stresstensor& Stre,
                                      const straintensor& Stra,
                                      const MaterialParameter& MaterialParameter_in) const;

        virtual const tensor& Dm_Dkin(const stresstensor& Stre,
                                      const straintensor& Stra,
                                      const MaterialParameter& MaterialParameter_in) const;

        virtual const tensor& Dm_Dkin2(const stresstensor& Stre,
                                       const straintensor& Stra,
                                       const MaterialParameter& MaterialParameter_in) const;

    protected:

        static tensor  PF_tensorR2;
        static tensor  PF_tensorR4;
};

#endif

