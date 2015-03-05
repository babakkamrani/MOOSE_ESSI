#ifndef ESSI_MATERIAL
#define ESSI_MATERIAL

#include "FiniteStrainMaterial.h"
#include "DruckerPragerIH.h"
#include "NewTemplate3Dep.h"
#include "Isotropic_Elastic.h"
#include "DP_YF.h"
#include "DP_PF.h"
#include "Linear_Eeq.h"

class ESSI_Material;

template<>
InputParameters validParams<ESSI_Material>();

class ESSI_Material : public FiniteStrainMaterial
{
public:
  ESSI_Material(const std::string & name, InputParameters parameters);
  
  //Accelerated DP - IH
//   DruckerPragerIH** Mat_ESSI;
  
  //New Template
  NewTemplate3Dep** Mat_ESSI;


protected:

  /// Modulus of Elasticity
  Real _E_Moose;
  double E_ESSI;

  /// Poisson Ratio
  Real _v_Moose;
  double v_ESSI;
  
  /// k
  Real _k_Moose;
  double k_ESSI;

  /// H
  Real _H_Moose;
  double H_ESSI;

  /// density
  Real _rho_Moose;
  double rho_ESSI;

  /// initial confinig stress
  Real _initialconfiningstress_Moose;
  double initialconfiningstress_ESSI;

  /// maximum number of iterations 
  Real _maximum_number_of_iterations_Moose;
  int maximum_number_of_iterations_ESSI;

  /// tolerance
  Real _tolerance_1_Moose;
  Real tolerance_1_ESSI;
  
  /// tolerance
  Real _tolerance_2_Moose; 
  Real tolerance_2_ESSI; 
  
  
  virtual void initQpStatefulProperties();
  virtual void computeQpStress() ;

  MaterialProperty<RankTwoTensor> & _strain_rate;
  MaterialProperty<RankTwoTensor> & _strain_increment;
  MaterialProperty<RankTwoTensor> & _total_strain_old;
  MaterialProperty<RankTwoTensor> & _elastic_strain_old;
  MaterialProperty<RankTwoTensor> & _stress_old;
  MaterialProperty<RankTwoTensor> & _rotation_increment;
  MaterialProperty<RankTwoTensor> & _dfgrd ;
  
  /// plastic strain
  MaterialProperty<RankTwoTensor> & _plastic_strain;

  /// Old value of plastic strain
  MaterialProperty<RankTwoTensor> & _plastic_strain_old;

  /// internal parameters
  MaterialProperty<std::vector<Real> > & _intnl;

  /// old values of internal parameters
  MaterialProperty<std::vector<Real> > & _intnl_old;

  /// yield functions
  MaterialProperty<std::vector<Real> > & _f;  
  
  
  

};

#endif //
