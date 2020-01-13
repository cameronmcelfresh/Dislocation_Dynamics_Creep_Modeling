#include <DefectiveCrystal.h>
#include <DislocationNetwork.h>
#include <Eigen>

#ifndef test_H_
#define test_H_

using namespace model;

template <typename DislocationNetworkType>
void printstress(double pos[3], DislocationNetworkType& DN)
{
	//VectorDim positionVec(pos[0],pos[1],pos[2]);
	Eigen::Matrix<double,3,1> positionVec(pos[0],pos[1],pos[2]);
	//MatrixDim stressMat(MatrixDim::Zero());
	Eigen::Matrix<double,3,3> stressMat(Eigen::Matrix<double,3,3>::Zero());
	
	//stress(positionVec) = stressMat; 

	stressMat+=DN->stress(positionVec);

	std::cout << stressMat << std::endl;
}

#endif
