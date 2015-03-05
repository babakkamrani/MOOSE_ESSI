#include "ESSI_Material.h"
#include <math.h> // for M_PI
#include "../../../include/utils/Moose_ESSI_Tensor_converter/MOOSE_ESSI_Tensor_Converter.h"
#include <DruckerPragerIH.h>
#include "MOOSE_ESSI_Tensor_Converter.h"
#include "Transient.h"
#include "NewTemplate3Dep.h"


extern Transient  *g_tensor_mechanics_executioner;


template<>
InputParameters validParams<ESSI_Material>()
{
  InputParameters params = validParams<FiniteStrainMaterial>();
  params.addParam<Real>("E_in", "Elastic Modulus");
  params.addParam<Real>("v_in", "Poisson Ratio ");
  params.addParam<Real>("k_in", "");
  params.addParam<Real>("H_in", "");
  params.addParam<Real>("rho_in", "");
  params.addParam<Real>("initialconfiningstress_in", "");
  params.addParam<Real>("maximum_number_of_iterations", "");
  params.addParam<Real>("tolerance_1", "");
  params.addParam<Real>("tolerance_2", "");  
  params.addClassDescription("Associated Drucker Prager with isothropic hardening");

  return params;
}

ESSI_Material::ESSI_Material(const std::string & name,
                                                         InputParameters parameters) :
    FiniteStrainMaterial(name, parameters),
    _E_Moose(getParam<Real>("E_in")),
    _v_Moose(getParam<Real>("v_in")),
    _k_Moose(getParam<Real>("k_in")),
    _H_Moose(getParam<Real>("H_in")),
    _rho_Moose(getParam<Real>("rho_in")),
    _initialconfiningstress_Moose(getParam<Real>("initialconfiningstress_in")),
    _maximum_number_of_iterations_Moose(getParam<Real>("maximum_number_of_iterations")),
    _tolerance_1_Moose(getParam<Real>("tolerance_1")),
    _tolerance_2_Moose(getParam<Real>("tolerance_2")),
    _strain_rate(declareProperty<RankTwoTensor>("strain_rate")),
    _strain_increment(declareProperty<RankTwoTensor>("strain_increment")),
    _total_strain_old(declarePropertyOld<RankTwoTensor>("total_strain")),
    _elastic_strain_old(declarePropertyOld<RankTwoTensor>("elastic_strain")),
    _stress_old(declarePropertyOld<RankTwoTensor>("stress")),
    _rotation_increment(declareProperty<RankTwoTensor>("rotation_increment")),
    _dfgrd(declareProperty<RankTwoTensor>("deformation gradient")),
    
    _plastic_strain(declareProperty<RankTwoTensor>("plastic_strain")),
    _plastic_strain_old(declarePropertyOld<RankTwoTensor>("plastic_strain")),
    _intnl(declareProperty<std::vector<Real> >("plastic_internal_parameter")),
    _intnl_old(declarePropertyOld<std::vector<Real> >("plastic_internal_parameter")),
    _f(declareProperty<std::vector<Real> >("plastic_yield_function"))    
    

{
  
    E_ESSI = (double) _E_Moose;
    v_ESSI = (double) _v_Moose;
    k_ESSI = (double) _k_Moose;
    H_ESSI = (double) _H_Moose;
    rho_ESSI = (double) _rho_Moose;
    initialconfiningstress_ESSI = (double) _initialconfiningstress_Moose;
    maximum_number_of_iterations_ESSI = (int) _maximum_number_of_iterations_Moose;
    tolerance_1_ESSI = (double) _tolerance_1_Moose;
    tolerance_2_ESSI = (double) _tolerance_2_Moose;
    


//Accelerated DP-IH
//-------------------------------------------------------------    
//     Mat_ESSI = new DruckerPragerIH * [8];
//     
//     for (int i = 0; i < 8; i++)
//     {
// 	Mat_ESSI[i] = new DruckerPragerIH(1, E_ESSI, v_ESSI, k_ESSI, H_ESSI, rho_ESSI, initialconfiningstress_ESSI, maximum_number_of_iterations_ESSI, tolerance_1_ESSI, tolerance_2_ESSI);
//     }
//-------------------------------------------------------------    

//NewTemplate3Dep
//-------------------------------------------------------------
    int Algorithm =  1; //Explicit (=0) or Implicit (=1)
    int number_of_subincrements = 10;
    
    
    Mat_ESSI = new NewTemplate3Dep * [8];
    
    for (int i = 0; i < 8; i++)
    {
      stresstensor zeroT;
      stresstensor initStress;
      initStress.val(1, 1) = -initialconfiningstress_ESSI;
      initStress.val(2, 2) = -initialconfiningstress_ESSI;
      initStress.val(3, 3) = -initialconfiningstress_ESSI;

      double MC[5] = {rho_ESSI, E_ESSI, v_ESSI, k_ESSI, H_ESSI};
      double IS[1] = {k_ESSI};
      MaterialParameter matpar(MC, 5, IS, 1);
      Isotropic_Elastic le(2, 3, initStress);
      DP_YF dpy(1, 1);
      DP_PF dpf(1, 1);
      Linear_Eeq Eep(5);
      ScalarEvolution* SE = {&Eep};

      Mat_ESSI[i] = new NewTemplate3Dep(1, &matpar, &le, &dpy, &dpf, &SE, Algorithm, number_of_subincrements, maximum_number_of_iterations_ESSI, tolerance_1_ESSI, tolerance_2_ESSI);
    }
//-------------------------------------------------------------
    
    
    

//   //Initialize the stress from confining stress ... :
//   //-------------------------------------------------- 
//    RankTwoTensor temp_stress((-1) * _initialconfiningstress_Moose, 
// 			     (-1) * _initialconfiningstress_Moose, 
// 			     (-1) * _initialconfiningstress_Moose,
// 			    0.0, 0.0, 0.0);
//   
//    _stress_old[_qp] = temp_stress;
//    _stress[_qp]     = temp_stress;
//  //-------------------------------------------------- 
    
    
}

void ESSI_Material::initQpStatefulProperties()
{
  
  FiniteStrainMaterial::initQpStatefulProperties();

  _plastic_strain[_qp].zero();
  _plastic_strain_old[_qp].zero();

  _intnl[_qp].assign(1, 0);
  _intnl_old[_qp].assign(1, 0);

  _f[_qp].assign(1, 0);
  
}


// void ESSI_Material::computeQpStress()
// {
// //   Update strain in intermediate configuration
//   _total_strain[_qp] = _total_strain_old[_qp] + _strain_increment[_qp];
// 
// //   Rotate strain to current configuration
//   _total_strain[_qp] = _rotation_increment[_qp] * _total_strain[_qp] * _rotation_increment[_qp].transpose();
// 
// //   For elastic problems elastic strain = total strain
//   _elastic_strain[_qp]=_total_strain[_qp];
// 
// //   stress = C * e
//   _stress[_qp] = _stress_old[_qp] + _elasticity_tensor[_qp] * _strain_increment[_qp]; //Calculate stress in intermediate configruation
// 
// //   Rotate the stress to the current configuration
//   _stress[_qp] = _rotation_increment[_qp] * _stress[_qp] * _rotation_increment[_qp].transpose();
// 
//   if (_qp == 0)
//   {
//    std::cout << "ESSI_Material::computeQpStress() AFTER---           _stress[" << _qp << "](0,0) =" << _stress[_qp](0,0) << "\n"; 
//    std::cout << "ESSI_Material::computeQpStress() AFTER---           _stress[" << _qp << "](2,2) =" << _stress[_qp](2,2) << "\n";
//    std::cout << "ESSI_Material::computeQpStress() AFTER--- _strain_increment[" << _qp << "] =" << _strain_increment[_qp](2,2) << "\n";
//    std::cout << "ESSI_Material::computeQpStress() AFTER---     _total_strain[" << _qp << "] =" << _total_strain[_qp](2,2) << "\n";
//    std::cout << "ESSI_Material::computeQpStress() AFTER--- _total_strain_old[" << _qp << "] =" << _total_strain_old[_qp](2,2) << "\n";
// 
// //    std::cout << "ESSI_Material::computeQpStress() AFTER--- stress_ESSI_temp.cval(3,3) = " <<stress_ESSI_temp.cval(3,3) << "\n";
//   }
// 
// 
// }




void ESSI_Material::computeQpStress()
{
   MOOSE_ESSI_Tensor_Converter Tensor_Converter; 
  
   straintensor strain_ESSI_temp_1;
   strain_ESSI_temp_1 = Tensor_Converter.RankTwoTensor_MOOSE_2_ESSI(_strain_increment[_qp]);  

   // Apply the prescribed strain to the material:
   Mat_ESSI[_qp]->setTrialStrainIncr(strain_ESSI_temp_1);   
   
   stresstensor stress_ESSI_temp;
   stress_ESSI_temp = Mat_ESSI[_qp]->getStressTensor();
   _stress[_qp] = Tensor_Converter.RankTwoTensor_ESSI_2_MOOSE(stress_ESSI_temp); 
   
   
//    std::cout << stress_ESSI_temp.cval(3,3)<< "\n";
//    std::cout << _stress[_qp](2,2)<< "\n";
   
   
   
   
// Push to the inside the cone ... small confining stress:
    //////////////////////////////////////////////////////////////////////
    _stress[_qp](0,0) = _stress[_qp](0,0) +_initialconfiningstress_Moose;
    _stress[_qp](1,1) = _stress[_qp](1,1) +_initialconfiningstress_Moose;
    _stress[_qp](2,2) = _stress[_qp](2,2) +_initialconfiningstress_Moose;
    //////////////////////////////////////////////////////////////////////
   
   straintensor strain_ESSI_temp_4;
   strain_ESSI_temp_4     = Mat_ESSI[_qp]->getPlasticStrainTensor();
   _plastic_strain[_qp] = Tensor_Converter.RankTwoTensor_ESSI_2_MOOSE(strain_ESSI_temp_4);
   
   
   

   
//    exit(1);
   

   
  //Update measures of strain
  ////////////////////////////////////////////////////////////////////// 
  _elastic_strain[_qp] = _elastic_strain_old[_qp] + _strain_increment[_qp] - (_plastic_strain[_qp] - _plastic_strain_old[_qp]);
  _total_strain[_qp]   = _total_strain_old[_qp]   + _strain_increment[_qp];
  //////////////////////////////////////////////////////////////////////
  
  
  //Rotate the tensors to the current configuration
  //////////////////////////////////////////////////////////////////////
  _stress[_qp]         = _rotation_increment[_qp] *_stress[_qp]          * _rotation_increment[_qp].transpose();
  _elastic_strain[_qp] = _rotation_increment[_qp] * _elastic_strain[_qp] * _rotation_increment[_qp].transpose();
  _total_strain[_qp]   = _rotation_increment[_qp] * _total_strain[_qp]   * _rotation_increment[_qp].transpose();
  _plastic_strain[_qp] = _rotation_increment[_qp] * _plastic_strain[_qp] * _rotation_increment[_qp].transpose();
  //////////////////////////////////////////////////////////////////////

//   if (_qp == 0)
//   {
//    std::cout << "ESSI_Material::computeQpStress() AFTER---           _stress[" << _qp << "](0,0) =" << _stress[_qp](0,0) << "\n"; 
//    std::cout << "ESSI_Material::computeQpStress() AFTER---           _stress[" << _qp << "](2,2) =" << _stress[_qp](2,2) << "\n";
//    std::cout << "ESSI_Material::computeQpStress() AFTER--- _strain_increment[" << _qp << "]      =" << _strain_increment[_qp](2,2) << "\n";
//    std::cout << "ESSI_Material::computeQpStress() AFTER---     _total_strain[" << _qp << "]      =" << _total_strain[_qp](2,2) << "\n";
//    std::cout << "ESSI_Material::computeQpStress() AFTER--- _total_strain_old[" << _qp << "]      =" << _total_strain_old[_qp](2,2) << "\n";
//   }

  if (_qp == 0)
  {
   std::cout << "ESSI_Material::computeQpStress() AFTER---           _stress " <<_total_strain[_qp](2,2) <<" " << _stress[_qp](2,2) << "\n"; 
  }
  

  //update old stress and strain:
  //////////////////////////////////////////////////////////////////////  
  _stress_old[_qp]         = _stress[_qp];
  _total_strain_old[_qp]   = _total_strain[_qp];
  _elastic_strain_old[_qp] = _elastic_strain[_qp];
  _plastic_strain_old[_qp] = _plastic_strain[_qp]; 
  ////////////////////////////////////////////////////////////////////// 
  
  
}

