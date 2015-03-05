#include "ESSI_Material.h"
#include <math.h> // for M_PI
#include "../../../include/utils/Moose_ESSI_Tensor_converter/MOOSE_ESSI_Tensor_Converter.h"
#include <DruckerPragerIH.h>
#include "MOOSE_ESSI_Tensor_Converter.h"
#include "Transient.h"

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
    
    
    
    Mat_ESSI = new DruckerPragerIH(1, E_ESSI, v_ESSI, k_ESSI, H_ESSI, rho_ESSI, initialconfiningstress_ESSI, 
		  maximum_number_of_iterations_ESSI, tolerance_1_ESSI, tolerance_2_ESSI);

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


void ESSI_Material::computeQpStress()
{
//   Update strain in intermediate configuration
  _total_strain[_qp] = _total_strain_old[_qp] + _strain_increment[_qp];

//   Rotate strain to current configuration
  _total_strain[_qp] = _rotation_increment[_qp] * _total_strain[_qp] * _rotation_increment[_qp].transpose();

//   For elastic problems elastic strain = total strain
  _elastic_strain[_qp]=_total_strain[_qp];

//   stress = C * e
  _stress[_qp] = _stress_old[_qp] + _elasticity_tensor[_qp] * _strain_increment[_qp]; //Calculate stress in intermediate configruation

//   Rotate the stress to the current configuration
  _stress[_qp] = _rotation_increment[_qp] * _stress[_qp] * _rotation_increment[_qp].transpose();

}




// void ESSI_Material::computeQpStress()
// {
//    MOOSE_ESSI_Tensor_Converter Tensor_Converter; 
//   
//    straintensor strain_ESSI_temp_1;
//    strain_ESSI_temp_1 = Tensor_Converter.RankTwoTensor_MOOSE_2_ESSI(_strain_increment[_qp]);  
// 
// //    Mat_ESSI->revertToLastCommit();   
// //   if (_qp == 0)
// //   {
// //    std::cout << "ESSI_Material::computeQpStress() BEFORE--- _stress[" << _qp << "] =" << _stress[_qp](2,2) << "\n";
// //    std::cout << "ESSI_Material::computeQpStress() BEFORE--- _strain_increment[" << _qp << "] =" << _strain_increment[_qp](2,2) << "\n";
// //    std::cout << "ESSI_Material::computeQpStress() BEFORE--- _total_strain[" << _qp << "] =" << _total_strain[_qp](2,2) << "\n";
// //   }
//    // Apply the prescribed strain to the material:
//    Mat_ESSI->setTrialStrainIncr(strain_ESSI_temp_1);   
// //    
// //    Mat_ESSI->commitState();
// //    Mat_ESSI->revertToStart();
// //    Mat_ESSI->   revertToLastCommit();  
//    
//    stresstensor stress_ESSI_temp;
//    stress_ESSI_temp = Mat_ESSI->getStressTensor();
//    _stress[_qp] = Tensor_Converter.RankTwoTensor_ESSI_2_MOOSE(stress_ESSI_temp);
//    
//    straintensor strain_ESSI_temp_2;
//    strain_ESSI_temp_2     = Mat_ESSI->getStrainTensor() - Mat_ESSI->getPlasticStrainTensor();
//    _elastic_strain[_qp]   = Tensor_Converter.RankTwoTensor_ESSI_2_MOOSE(strain_ESSI_temp_2);
//    
//    straintensor strain_ESSI_temp_3;
//    strain_ESSI_temp_3     = Mat_ESSI->getStrainTensor();
//    _total_strain[_qp]   = Tensor_Converter.RankTwoTensor_ESSI_2_MOOSE(strain_ESSI_temp_3);
//    
//    straintensor strain_ESSI_temp_4;
//    strain_ESSI_temp_4     = Mat_ESSI->getPlasticStrainTensor();
//    _plastic_strain[_qp] = Tensor_Converter.RankTwoTensor_ESSI_2_MOOSE(strain_ESSI_temp_4);
// // 
// //    Mat_ESSI->commitState();
// //    Mat_ESSI->   revertToLastCommit();   
// // _problem;
//    
// //    Mat_ESSI_old = Mat_ESSI->getCopy();
//    
//   //Update measures of strain
//   _elastic_strain[_qp] = _elastic_strain_old[_qp] + _strain_increment[_qp] - (_plastic_strain[_qp] - _plastic_strain_old[_qp]);
//   _total_strain[_qp] = _total_strain_old[_qp] + _strain_increment[_qp];
// 
//   //Rotate the tensors to the current configuration
//   _stress[_qp] = _rotation_increment[_qp]*_stress[_qp]*_rotation_increment[_qp].transpose();
//   _elastic_strain[_qp] = _rotation_increment[_qp] * _elastic_strain[_qp] * _rotation_increment[_qp].transpose();
//   _total_strain[_qp] = _rotation_increment[_qp] * _total_strain[_qp] * _rotation_increment[_qp].transpose();
//   _plastic_strain[_qp] = _rotation_increment[_qp] * _plastic_strain[_qp] * _rotation_increment[_qp].transpose();
//   
// 
//   if (_qp == 0)
//   {
//    std::cout << "ESSI_Material::computeQpStress() AFTER---           _stress[" << _qp << "] =" << _stress[_qp](2,2) << "\n";
//    std::cout << "ESSI_Material::computeQpStress() AFTER--- _strain_increment[" << _qp << "] =" << _strain_increment[_qp](2,2) << "\n";
//    std::cout << "ESSI_Material::computeQpStress() AFTER---     _total_strain[" << _qp << "] =" << _total_strain[_qp](2,2) << "\n";
//    std::cout << "ESSI_Material::computeQpStress() AFTER--- _total_strain_old[" << _qp << "] =" << _total_strain_old[_qp](2,2) << "\n";
// 
// //    std::cout << "ESSI_Material::computeQpStress() AFTER--- stress_ESSI_temp.cval(3,3) = " <<stress_ESSI_temp.cval(3,3) << "\n";
//   }
// 
//   //update old stress and strain:
//   //-------------------------------------------------  
//   _stress_old[_qp]         = _stress[_qp];
//   _total_strain_old[_qp]   = _total_strain[_qp];
//   _elastic_strain_old[_qp] = _elastic_strain[_qp];
// //   _plastic_strain_old[_qp] = _plastic_strain[_qp]; 
//   //-------------------------------------------------  
//   
//   
//    Mat_ESSI->revertToStart();
//   
//   
// }

