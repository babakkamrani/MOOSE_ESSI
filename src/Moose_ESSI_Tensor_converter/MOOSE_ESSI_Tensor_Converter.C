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



#include "MOOSE_ESSI_Tensor_Converter.h"
#include "BJtensor.h"
#include "stresst.h"


MOOSE_ESSI_Tensor_Converter::MOOSE_ESSI_Tensor_Converter()
{
}

// stresstensor
// MOOSE_ESSI_Tensor_Converter::RankTwoTensor_MOOSE_2_ESSI(RankTwoTensor MOOSE_Tensor)
// {
//   stresstensor ESSI_Tensor;
//   ESSI_Tensor.val(0,0) = MOOSE_Tensor(0,0);
//   ESSI_Tensor.val(1,0) = MOOSE_Tensor(1,0);
//   ESSI_Tensor.val(2,0) = MOOSE_Tensor(2,0);
//   
//   ESSI_Tensor.val(0,1) = MOOSE_Tensor(0,1);
//   ESSI_Tensor.val(1,1) = MOOSE_Tensor(1,1);
//   ESSI_Tensor.val(2,1) = MOOSE_Tensor(2,1);
//   
//   ESSI_Tensor.val(0,2) = MOOSE_Tensor(0,2);
//   ESSI_Tensor.val(1,2) = MOOSE_Tensor(1,2);
//   ESSI_Tensor.val(2,2) = MOOSE_Tensor(2,2);
//   
//   
//   std::cerr << "MOOSE_ESSI_Tensor_Converter::RankTwoTensor_MOOSE_2_ESSI :" <<  ESSI_Tensor.cval(0,0) << "\n";
// /*  std::cerr << "MOOSE_ESSI_Tensor_Converter::RankTwoTensor_MOOSE_2_ESSI :" <<  ESSI_Tensor.cval(1,0) << "\n";
//   std::cerr << "MOOSE_ESSI_Tensor_Converter::RankTwoTensor_MOOSE_2_ESSI :" <<  ESSI_Tensor.cval(2,0) << "\n";
//   
//   std::cerr << "MOOSE_ESSI_Tensor_Converter::RankTwoTensor_MOOSE_2_ESSI :" <<  ESSI_Tensor.cval(0,1) << "\n";
//   std::cerr << "MOOSE_ESSI_Tensor_Converter::RankTwoTensor_MOOSE_2_ESSI :" <<  ESSI_Tensor.cval(1,1) << "\n";
//   std::cerr << "MOOSE_ESSI_Tensor_Converter::RankTwoTensor_MOOSE_2_ESSI :" <<  ESSI_Tensor.cval(2,1) << "\n";
//   
//   std::cerr << "MOOSE_ESSI_Tensor_Converter::RankTwoTensor_MOOSE_2_ESSI :" <<  ESSI_Tensor.cval(0,2) << "\n";
//   std::cerr << "MOOSE_ESSI_Tensor_Converter::RankTwoTensor_MOOSE_2_ESSI :" <<  ESSI_Tensor.cval(1,2) << "\n"; 
//   std::cerr << "MOOSE_ESSI_Tensor_Converter::RankTwoTensor_MOOSE_2_ESSI :" <<  ESSI_Tensor.cval(2,3) << "\n"; */ 
//   
//   return ESSI_Tensor;
// }

BJtensor
MOOSE_ESSI_Tensor_Converter::RankTwoTensor_MOOSE_2_ESSI(RankTwoTensor MOOSE_Tensor)
{
  
//   BJtensor ESSI_Tensor(2, def_dim_2, 0.0);
 
  stresstensor ESSI_Tensor; 
//   double temp = MOOSE_Tensor(1,1);
//   std::cerr << "MOOSE_Tensor(0,0)" << temp <<"\n";
  
  
  ESSI_Tensor.val(1,1) = (double)MOOSE_Tensor(0,0);
  ESSI_Tensor.val(2,1) = (double)MOOSE_Tensor(1,0);
  ESSI_Tensor.val(3,1) = (double)MOOSE_Tensor(2,0);
  
  ESSI_Tensor.val(1,2) = (double)MOOSE_Tensor(0,1);
  ESSI_Tensor.val(2,2) = (double)MOOSE_Tensor(1,1);
  ESSI_Tensor.val(3,2) = (double)MOOSE_Tensor(2,1);
  
  ESSI_Tensor.val(1,3) = (double)MOOSE_Tensor(0,2);
  ESSI_Tensor.val(2,3) = (double)MOOSE_Tensor(1,2);
  ESSI_Tensor.val(3,3) = (double)MOOSE_Tensor(2,2);
  return ESSI_Tensor;
}




// RankTwoTensor
// MOOSE_ESSI_Tensor_Converter::RankTwoTensor_ESSI_2_MOOSE(stresstensor ESSI_Tensor)
// {
//   RankTwoTensor MOOSE_Tensor;
//   MOOSE_Tensor(0,0) = ESSI_Tensor.cval(0,0);
//   MOOSE_Tensor(1,0) = ESSI_Tensor.cval(1,0);
//   MOOSE_Tensor(2,0) = ESSI_Tensor.cval(2,0);
//   
//   MOOSE_Tensor(0,1) = ESSI_Tensor.cval(0,1);
//   MOOSE_Tensor(1,1) = ESSI_Tensor.cval(1,1);
//   MOOSE_Tensor(2,1) = ESSI_Tensor.cval(2,1);
//   
//   MOOSE_Tensor(0,2) = ESSI_Tensor.cval(0,2);
//   MOOSE_Tensor(1,2) = ESSI_Tensor.cval(1,2);
//   MOOSE_Tensor(2,2) = ESSI_Tensor.cval(2,2);
//   
//   return MOOSE_Tensor;
// }

RankTwoTensor
MOOSE_ESSI_Tensor_Converter::RankTwoTensor_ESSI_2_MOOSE(BJtensor ESSI_Tensor)
{
  RankTwoTensor MOOSE_Tensor;
  MOOSE_Tensor(0,0) = (Real)ESSI_Tensor.cval(1,1);
  MOOSE_Tensor(1,0) = (Real)ESSI_Tensor.cval(2,1);
  MOOSE_Tensor(2,0) = (Real)ESSI_Tensor.cval(3,1);
  
  MOOSE_Tensor(0,1) = (Real)ESSI_Tensor.cval(1,2);
  MOOSE_Tensor(1,1) = (Real)ESSI_Tensor.cval(2,2);
  MOOSE_Tensor(2,1) = (Real)ESSI_Tensor.cval(3,2);
  
  MOOSE_Tensor(0,2) = (Real)ESSI_Tensor.cval(1,3);
  MOOSE_Tensor(1,2) = (Real)ESSI_Tensor.cval(2,3);
  MOOSE_Tensor(2,2) = (Real)ESSI_Tensor.cval(3,3);
  
  return MOOSE_Tensor;
}



BJtensor
MOOSE_ESSI_Tensor_Converter::RankFourTensor_MOOSE_2_ESSI(RankFourTensor MOOSE_Tensor)
{
  BJtensor ESSI_Tensor;
  for(unsigned int i = 0; i < 3; i++)
  {
    for(unsigned int j = 0; j < 3; j++)
    {
      for(unsigned int k = 0; k < 3; k++)
      {
	for(unsigned int l = 0; l < 3; l++)
	{
	  ESSI_Tensor.val(i,j,k,l) = MOOSE_Tensor(i,j,k,l);
	}
      }
    }
  }
  
  return ESSI_Tensor;
}


RankFourTensor
MOOSE_ESSI_Tensor_Converter::RankFourTensor_ESSI_2_MOOSE( BJtensor ESSI_Tensor)
{
  RankFourTensor MOOSE_Tensor;
  for(unsigned int i = 0; i < 3; i++)
  {
    for(unsigned int j = 0; j < 3; j++)
    {
      for(unsigned int k = 0; k < 3; k++)
      {
	for(unsigned int l = 0; l < 3; l++)
	{
	  MOOSE_Tensor(i,j,k,l) = ESSI_Tensor.cval(i,j,k,l);
	}
      }
    }
  }
  return MOOSE_Tensor;
	  
}




