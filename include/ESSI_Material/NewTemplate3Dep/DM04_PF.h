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
// PROGRAMMER:        Zhao Cheng
// Note:              Helpful discuss with Mahdi Taiebat and Professor Y.F. Dafalias
// DATE:              Fall 2005
// UPDATE HISTORY:    Guanzhou Jie updated for parallel, Dec. 2006
//                    Nima Tafazzoli updated for API (Feb 2009)
//
///////////////////////////////////////////////////////////////////////////////
//

#ifndef DM04_PF_H
#define DM04_PF_H

#include "PlasticFlow.h"
#include <math.h>
#include "classTags.h"

class DM04_PF : public PlasticFlow
{
    public:

        //! rf: Dafalias & Manzari 2004
        //! Plastic Flow function for Dafalias-Manzari:
        //! inputs:
        //! - e0_which: to define the material parameter type (initial void ratio, =0)
        //! - index_e0: to locate the position in the defined (constant) MaterialParameter command for e0
        //! - e_r_which_in: to define the material parameter type (e_r, =0)
        //! - index_e_r_in: to locate the position in the defined (constant) MaterialParameter command for e_r
        //! - lambda_c_which_in: to define the material parameter type (c, =0)
        //! - index_lambda_c_in: to locate the position in the defined (constant) MaterialParameter command for c
        //! - xi_which_in: to define the material parameter type (xi, =0)
        //! - index_xi_in: to locate the position in the defined (constant) MaterialParameter command for xi
        //! - Pat_which_in: to define the material parameter type (atmospheric pressure, =0)
        //! - index_Pat_in: to locate the position in the defined (constant) MaterialParameter command for Pat
        //! - m_which_in: to define the material parameter type (m, =0)
        //! - index_m_in: to locate the position in the defined (constant) MaterialParameter command for m
        //! - M_cal_which_in: to define the material parameter type (M, critial state line slope, =0)
        //! - index_M_cal_in: to locate the position in the defined (constant) MaterialParameter command for M
        //! - cc_which_in: to define the material parameter type (cc, =0)
        //! - index_cc_in: to locate the position in the defined (constant) MaterialParameter command for cc
        //! - A0_which_in: to define the material parameter type (A0, =0)
        //! - index_A0_in: to locate the position in the defined (constant) MaterialParameter command for A0
        //! - nd_which_in: to define the material parameter type (nd, =0)
        //! - index_nd_in: to locate the position in the defined (constant) MaterialParameter command for nd
        //! - alpha_which_in: to define the material parameter type (back stress, =2)
        //! - index_alpha_in: to locate the position in the defined (constant) MaterialParameter command for alpha (back stress)
        //! - z_which_in: to define the material parameter type (fabric tensor, =2)
        //! - index_z_in: to locate the position in the defined (constant) MaterialParameter command for z (fabric tensor)
        //!

        DM04_PF() : PlasticFlow(PF_TAG_DM04) {}; //Guanzhou added for parallel processing

        DM04_PF(int e0_which, int index_e0,
                int e_r_which_in, int index_e_r_in,
                int lambda_c_which_in, int index_lambda_c_in,
                int xi_which_in, int index_xi_in,
                int Pat_which_in, int index_Pat_in,
                int m_which_in, int index_m_in,
                int M_cal_which_in, int index_M_cal_in,
                int cc_which_in, int index_cc_in,
                int A0_which_in, int index_A0_in,
                int nd_which_in, int index_nd_in,
                int alpha_which_in, int index_alpha_in,
                int z_which_in, int index_z_in);

        //! Deconstructor for Dafalias-Manzari plastic flow function
        ~DM04_PF();
        PlasticFlow* newObj();

        const char* getPlasticFlowType(void) const
        {
            return "DM04_PF";
        };


        //! PlasticFlowTensor: to define the plastic flow
        //! inputs:
        //! - stresstensor &Stre: the tensor of stress
        //! - straintensor &Stra: the tensor of strain
        //! - MaterialParameter &MaterialParameter_in: the class of material parameter

        const straintensor& PlasticFlowTensor(const stresstensor& Stre,
                                              const straintensor& Stra,
                                              const MaterialParameter& MaterialParameter_in) const;

        //! Dm_Ds: to obtain the derivative of the plastic flow respect to the scalar internal varible
        //! inputs:
        //! - stresstensor &Stre: the tensor of stress
        //! - straintensor &Stra: the tensor of strain
        //! - MaterialParameter &MaterialParameter_in: the class of material parameter

        const tensor& Dm_Ds(const stresstensor& Stre,
                            const straintensor& Stra,
                            const MaterialParameter& MaterialParameter_in) const;

        //! Dm_Dkin: to obtain the derivative of the plastic flow respect to the 1st tensor internal varible (back stress)
        //! inputs:
        //! - stresstensor &Stre: the tensor of stress
        //! - straintensor &Stra: the tensor of strain
        //! - MaterialParameter &MaterialParameter_in: the class of material parameter

        const tensor& Dm_Dkin(const stresstensor& Stre,
                              const straintensor& Stra,
                              const MaterialParameter& MaterialParameter_in) const;

        //! Dm_Dkin2: to obtain the derivative of the plastic flow respect to the 2nd tensor internal varible (fabric tensor)
        //! inputs:
        //! - stresstensor &Stre: the tensor of stress
        //! - straintensor &Stra: the tensor of strain
        //! - MaterialParameter &MaterialParameter_in: the class of material parameter

        const tensor& Dm_Dkin2(const stresstensor& Stre,
                               const straintensor& Stra,
                               const MaterialParameter& MaterialParameter_in) const;


    private:

        //! rf: Dafalias & Manzari 2004
        //! gete0: to get e0
        //! inputs:
        //! - MaterialParameter &MaterialParameter_in: the class of material parameter
        double gete0(const MaterialParameter& MaterialParameter_in) const;

        //! rf: Dafalias & Manzari 2004
        //! gete_r: to get e_r
        //! inputs:
        //! - MaterialParameter &MaterialParameter_in: the class of material parameter
        double gete_r(const MaterialParameter& MaterialParameter_in) const;

        //! rf: Dafalias & Manzari 2004
        //! getlambda_c: to get lambda_c
        //! inputs:
        //! - MaterialParameter &MaterialParameter_in: the class of material parameter
        double getlambda_c(const MaterialParameter& MaterialParameter_in) const;

        //! rf: Dafalias & Manzari 2004
        //! getxi: to get xi
        //! inputs:
        //! - MaterialParameter &MaterialParameter_in: the class of material parameter
        double getxi(const MaterialParameter& MaterialParameter_in) const;

        //! rf: Dafalias & Manzari 2004
        //! getPat: to get pat (constant of atmosphereic pressure)
        //! inputs:
        //! - MaterialParameter &MaterialParameter_in: the class of material parameter
        double getPat(const MaterialParameter& MaterialParameter_in) const;

        //! rf: Dafalias & Manzari 2004
        //! getm: to get m (constant in yield function)
        //! inputs:
        //! - MaterialParameter &MaterialParameter_in: the class of material parameter
        double getm(const MaterialParameter& MaterialParameter_in) const;

        //! rf: Dafalias & Manzari 2004
        //! getM_cal: to get M_cal (critial state line slope)
        //! inputs:
        //! - MaterialParameter &MaterialParameter_in: the class of material parameter
        double getM_cal(const MaterialParameter& MaterialParameter_in) const;

        //! rf: Dafalias & Manzari 2004
        //! getcc: to get cc
        //! inputs:
        //! - MaterialParameter &MaterialParameter_in: the class of material parameter
        double getcc(const MaterialParameter& MaterialParameter_in) const;

        //! rf: Dafalias & Manzari 2004
        //! getA0: to get A0
        //! inputs:
        //! - MaterialParameter &MaterialParameter_in: the class of material parameter
        double getA0(const MaterialParameter& MaterialParameter_in) const;

        //! rf: Dafalias & Manzari 2004
        //! getnd: to get nd
        //! inputs:
        //! - MaterialParameter &MaterialParameter_in: the class of material parameter
        double getnd(const MaterialParameter& MaterialParameter_in) const;

        //! rf: Dafalias & Manzari 2004
        //! getalpha: to get alpha
        //! inputs:
        //! - MaterialParameter &MaterialParameter_in: the class of material parameter
        const stresstensor& getalpha(const MaterialParameter& MaterialParameter_in) const;

        //! rf: Dafalias & Manzari 2004
        //! getz: to get z
        //! inputs:
        //! - MaterialParameter &MaterialParameter_in: the class of material parameter
        const stresstensor& getz(const MaterialParameter& MaterialParameter_in) const;

        //! rf: Dafalias & Manzari 2004
        //! getParameters:
        //! inputs:
        //! - MaterialParameter &MaterialParameter_in: the class of material parameter
        //! - parIndex_in, which_in

        inline double getParameters(const MaterialParameter& MaterialParameter_in, int parIndex_in, int which_in) const;

        //! rf: Dafalias & Manzari 2004
        //! getec: to get e_c
        //! inputs:
        //! - e_r, lambda_c, xi, Pat,p_c
        inline double getec(double e_r, double lambda_c, double xi, double Pat, double p_c) const;

        //! rf: Dafalias & Manzari 2004
        //! getg: to get g(theta)
        //! inputs:
        //! - c , cos3theta
        inline double getg(double c, double cos3theta) const;


    private:

        int e0_which;
        int index_e0;
        int e_r_which;
        int index_e_r;
        int lambda_c_which;
        int index_lambda_c;
        int xi_which;
        int index_xi;
        int Pat_which;
        int index_Pat;
        int m_which;
        int index_m;
        int M_cal_which;
        int index_M_cal;
        int cc_which;
        int index_cc;
        int A0_which;
        int index_A0;
        int nd_which;
        int index_nd;

        int alpha_which;
        int index_alpha;
        int z_which;
        int index_z;

        static straintensor DM04m;
        static stresstensor DM04temp;
};


#endif


