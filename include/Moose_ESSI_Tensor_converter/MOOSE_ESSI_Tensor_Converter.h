//-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
// Programmer : Babak Kamrani (babakkamtrani@yahoo.com)                             -/
//-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
// Date       : November 15, 2014                                                   -/
// Purpose : Converts the tensor of ESSI to Moose's tensor mechanics format         -/
//-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
//-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
//-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
//-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
//-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
//-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
//-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/


#ifndef MOOSE_ESSI_Tensor_Converter_H
#define MOOSE_ESSI_Tensor_Converter_H

#include "BJtensor.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "stresst.h"



class MOOSE_ESSI_Tensor_Converter
{
public:
  
  MOOSE_ESSI_Tensor_Converter();
  
  // reads Moose rank two tensor and fills the BJtensor(BJtensor) one accordingly
//   stresstensor RankTwoTensor_MOOSE_2_ESSI(RankTwoTensor  MOOSE_Tensor);
  BJtensor RankTwoTensor_MOOSE_2_ESSI(RankTwoTensor  MOOSE_Tensor);
  
  // reads ESSI's(BJtensor) rank two tensor and fills the Moose one accordingly
//   RankTwoTensor RankTwoTensor_ESSI_2_MOOSE(stresstensor ESSI_Tensor);
  RankTwoTensor RankTwoTensor_ESSI_2_MOOSE(BJtensor ESSI_Tensor);
  
  // reads Moose rank four tensor and fills the BJtensor(BJtensor) accordingly
  BJtensor RankFourTensor_MOOSE_2_ESSI(RankFourTensor MOOSE_Tensor);

  // reads ESSI's(BJtensor) rank four tensor and fills the Moose accordingly
  RankFourTensor RankFourTensor_ESSI_2_MOOSE( BJtensor ESSI_Tensor);
  
  
//   stresstensor * ESSI_Stress_Tensor;

};
#endif




