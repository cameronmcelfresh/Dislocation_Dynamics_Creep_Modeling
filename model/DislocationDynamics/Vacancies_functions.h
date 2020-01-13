#include <iostream> 
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <complex>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <DefectiveCrystal.h>
#include <DislocationNetwork.h>
#include <Vacancy_Inputs.h>
#include <string>
#include <cstdio>

#ifndef vacancy_functions_H_
#define vacancy_functions_H_

using namespace std;
using namespace Eigen;


void princStress(double stressTens[][3], double princpStresses[])
//Function to take a given stress tensor and calculate the principle stresses and eigenvalues
{
	MatrixXd m(3,3);
 	m(0,0) = stressTens[0][0];
 	m(0,1) = stressTens[0][1];
 	m(0,2) = stressTens[0][2];
 	m(1,0) = m(0,1);
 	m(1,1) = stressTens[1][1];
 	m(1,2) = stressTens[1][2];
 	m(2,0) = m(0,2);
 	m(2,1) = m(1,2);
 	m(2,2) = stressTens[2][2];

	EigenSolver<MatrixXd> es(m);

	Vector3cd tempVec =  es.eigenvalues();

	princpStresses[0] = tempVec(0).real();
	princpStresses[1] = tempVec(1).real();
	princpStresses[2] = tempVec(2).real();
}

double determinent(double x1, double x2, double y1, double y2)
//Function to calculate the determinent of a 2x2 matrix, given the values at their respective coordinates
{
	return x1*y2-x2*y1;
}


void copy3x3(double array2Copy[][3], double destArray[][3])
//Function to copy a 3x3 array to another array
{
	for(int i = 0; i<3; i++)
	{
		for(int j=0;j<3; j++)
		{
			destArray[i][j] = array2Copy[i][j];
		}
	}
}

void princDir(double appliedStress[][3], double princpStresses[], double dirCosinesMat[][3])
//Function to calculate the principle directions on a point given the applied stress tensor and principle stresses
//Follows https://www.eoas.ubc.ca/courses/eosc433/lecture-material/StressStrain-Review.pdf
{
	MatrixXd m(3,3);
 	m(0,0) = appliedStress[0][0];
 	m(0,1) = appliedStress[0][1];
 	m(0,2) = appliedStress[0][2];
 	m(1,0) = m(0,1);
 	m(1,1) = appliedStress[1][1];
 	m(1,2) = appliedStress[1][2];
 	m(2,0) = m(0,2);
 	m(2,1) = m(1,2);
 	m(2,2) = appliedStress[2][2];

	EigenSolver<MatrixXd> es(m);

	MatrixXd tempEigenMat = es.eigenvectors().real();

	for(int i = 0; i<3; i++)
		for(int j = 0; j<3; j++)
			dirCosinesMat[i][j] = tempEigenMat(i,j);

}


void printmatrix(double A[][4]) 
//Function to print a 3x3 array
{
    int p=3;
    int q=4;

    for (int i=0; i<p; i++) {
            for (int j=0; j<q; j++) {
                    std::cout << setw(7) << setprecision(4) << A[i][j] << " ";
            }
            std::cout << endl;
    }

    std::cout << endl;
}

void RowReduce(double A[][4])
//Function to row reduce a 3x4 matrix to row-reduced echelon form
{
    const int nrows = 3; // number of rows
    const int ncols = 4; // number of columns

    int lead = 0; 

    while (lead < nrows) 
    {
        double d, m;

        for (int r = 0; r < nrows; r++) 

        	{ // for each row ...
            /* calculate divisor and multiplier */
            d = A[lead][lead];
            m = A[r][lead] / A[lead][lead];

            for (int c = 0; c < ncols; c++) 
            { // for each column ...
                if (r == lead)
                    A[r][c] /= d;               // make pivot = 1
                else
                    A[r][c] -= A[lead][c] * m;  // make other = 0
            }
        }

        lead++;
        printmatrix(A);
}
}

double thetaCalc(double vec1[], double vec2[])
//Function to calculate the angle between 2 vectors
{
	double dotProduct = vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2];

	double sizeVec1 = sqrt(pow(vec1[0],2) + pow(vec1[1],2) + pow(vec1[2],2));
	double sizeVec2 = sqrt(pow(vec2[0],2) + pow(vec2[1],2) + pow(vec2[2],2));

	return acos(dotProduct/(sizeVec1*sizeVec2));
}


void print3x3(double dirCosines[][3])
//Function to print direction cosines 3x3 matrix (or any 3x3 array matrix)
{
	for(int i = 0; i<3; i++)
	{
		for(int j=0;j<3; j++)
		{
			std::cout << dirCosines[i][j] << "\t";
		}
		std::cout << endl;
	}
}

void rotateToNewAxis(double dirCosinesMat[][3], double origVec[])
//Funtion to adjust a vector a new 3D basis set, dictated by the transformation matrix
{
	
	/*
	cout << "Original Vec: ";
	for(int i = 0; i<3; i++)
		cout << origVec[i] << " ";
	cout << endl;
	*/

	double tempVec[3];
	tempVec[0] = origVec[0];
	tempVec[1] = origVec[1];
	tempVec[2] = origVec[2];

	for(int i=0; i<3;i++)
		origVec[i] = tempVec[0]*dirCosinesMat[i][0]+ tempVec[1]*dirCosinesMat[i][1]+ tempVec[2]*dirCosinesMat[i][2];

	/*cout << "New Vec: ";
	for(int i = 0; i<3; i++)
		cout << origVec[i] << " ";
	cout << endl;
	*/
}

void rotateBackToOld(double dirCosinesMat[][3], double origVec[])
//Funtion to adjust a vector a new 3D basis set, dictated by the transformation matrix
{

	double tempVec[3];
	tempVec[0] = origVec[0];
	tempVec[1] = origVec[1];
	tempVec[2] = origVec[2];

	double dirCosinesTemp[3][3];

	copy3x3(dirCosinesMat, dirCosinesTemp); //Copy Matrix

	//print3x3(dirCosinesTemp);

	//Transpose Matix!
	dirCosinesTemp[0][1] = dirCosinesMat[1][0];
	dirCosinesTemp[1][0] = dirCosinesMat[0][1];
	dirCosinesTemp[0][2] = dirCosinesMat[2][0];
	dirCosinesTemp[2][0] = dirCosinesMat[0][2];
	dirCosinesTemp[1][2] = dirCosinesMat[2][1];
	dirCosinesTemp[2][1] = dirCosinesMat[1][2];
	//print3x3(dirCosinesTemp);


	for(int i=0; i<3;i++)
		origVec[i] = tempVec[0]*dirCosinesTemp[i][0]+ tempVec[1]*dirCosinesTemp[i][1]+ tempVec[2]*dirCosinesTemp[i][2];
}



void positionToAppliedStress(double position[], double appliedStress[][3])
//Function to return the stress tensor at a given point some x,y distance from an edge dislocation
//Consider edge dislocation existing at point 0,0
{
	double x = position[0];
	double y = position[1];

	double stressConstant = G*b/(2*3.14159*(1-poissonRatio));

	appliedStress[0][0] = -stressConstant*y*(3*x*x + y*y)/pow((x*x + y*y),2);
	appliedStress[1][1] = stressConstant*y*(x*x - y*y)/pow((x*x + y*y),2);
	appliedStress[2][2] = -stressConstant*2*y/(x*x + y*y);


	appliedStress[0][1] = stressConstant*x*(x*x - y*y)/pow((x*x + y*y),2);
	appliedStress[1][0] = stressConstant*x*(x*x - y*y)/pow((x*x + y*y),2);

	appliedStress[0][2] = 0;
	appliedStress[1][2] = 0;
	appliedStress[2][0] = 0;
	appliedStress[2][1] = 0;
}


template <typename DislocationNetworkType>
void DDpositionToAppliedStress(double pos[3], DislocationNetworkType& DN,double appliedStress[][3])
{
	//VectorDim positionVec(pos[0],pos[1],pos[2]);
	Eigen::Matrix<double,3,1> positionVec(pos[0]/b,pos[1]/b,pos[2]/b);
	//MatrixDim stressMat(MatrixDim::Zero());
	Eigen::Matrix<double,3,3> stressMat(Eigen::Matrix<double,3,3>::Zero());
	
	//stress(positionVec) = stressMat; 
	//std::cout<< "Getting Stress from DD..." << std::endl;
	stressMat+=DN->stress(positionVec);
	//std::cout<< "Stress retrieved from DD..." << std::endl;
	//std::cout << stressMat << std::endl;

	for(int i=0;i<3;i++)
		for(int j = 0; j<3; j++)
			appliedStress[i][j]=stressMat(i,j)*G;
}

//template <typename DislocationStructure>
void straightArmStress(Eigen::Matrix<double,3,1> pStart, Eigen::Matrix<double,3,1> pEnd,Eigen::Matrix<double,3,1> burgers,Eigen::Matrix<double,3,1> position, Eigen::Matrix<double,3,3>& tempStress)
{
	using namespace model;

	StressStraight straightSegment = StressStraight(pStart, pEnd, burgers); //Create a dislocation arm and calculate the stress from it given a position

	tempStress = straightSegment.stress(position);

}

template <typename DislocationStructure>
void DSpositionToAppliedStress(double pos[3], DislocationStructure DS,double appliedStress[][3])
{
	using namespace model;


	Eigen::Matrix<double,3,1> positionVec(pos[0]/b,pos[1]/b,pos[2]/b);

	Eigen::Matrix<double,3,3> stressMat(Eigen::Matrix<double,3,3>::Zero());

	Eigen::Matrix<double,3,3> tempMat;
	Eigen::Matrix<double,3,1> p0;
	Eigen::Matrix<double,3,1> p1;
	Eigen::Matrix<double,3,1> burgers;



/*
	//Attempt #1
	for( auto const& [key, val] : DS.nodeConnectivity )
	{

	//Ignore the forces of all connecting segments to simulate an "infinite" dislocation
	//if (val[1]==1)
	//	continue;


	//Node 0 Positions
	p0(0) = DS.nodePositions[key[0]][0];
	p0(1) = DS.nodePositions[key[0]][1];
	p0(2) = DS.nodePositions[key[0]][2];

	//std::cout << "p0 " << p0 << std::endl;
	
	//Node 1 Positions
	p1(0) = DS.nodePositions[key[1]][0];
	p1(1) = DS.nodePositions[key[1]][1];
	p1(2) = DS.nodePositions[key[1]][2];

	//std::cout << "p1 " << p1 << std::endl;

	burgers(0) = DS.loops[val[0]][0];
	burgers(1) = DS.loops[val[0]][1];
	burgers(2) = DS.loops[val[0]][2];

	//std::cout << "burgers " << std::endl;
	//std::cout << burgers << std::endl;

	std::cout << "Calculating stress contribution from node " << key[0] << " ("<< p0(0) <<"," << p0(1) << "," << p0(2) << ") to " << key[1] << " ("<< p1(0) <<"," << p1(1) << "," << p1(2) << ")" << std::endl;

	//straightArmStress(p0, p1, burgers, positionVec, tempMat); //Call a stress function to return the stress give a dislocation start/ending point and a position point


	Eigen::Matrix<double,3,1> p2;
	//Node 3 Positions
	p2(0) = DS.nodePositions[key[1]+1][0];
	p2(1) = DS.nodePositions[key[1]+1][1];
	p2(2) = DS.nodePositions[key[1]+1][2];
	StressStraight straightSegment1 = StressStraight(p1, p2, burgers); //Create a dislocation arm and calculate the stress from it given a position
	std::cout << "Tempmat +1" << std::endl;
	Eigen::Matrix<double,3,3> tempMat1;
	tempMat1 = straightSegment1.stress(positionVec);
	std::cout << tempMat1 << std::endl;

	///////////TEST////////////


	StressStraight straightSegment = StressStraight(p0, p1, burgers); //Create a dislocation arm and calculate the stress from it given a position


	//tempMat = straightSegment.stress(positionVec);
	tempMat = straightSegment.stress(positionVec);
	////////////////////////

	std::cout << "Tempmat" << std::endl;
	std::cout << tempMat << std::endl;

	//stressMat+=tempMat;

	stressMat=stressMat + tempMat;

	std::cout << "Cumulative stressmat" << std::endl;
	std::cout << stressMat << std::endl;

	//for(int i=0;i<3;i++)
	//	for(int j = 0; j<3; j++)
	//		tempMat(i,j) = 0;

	}
*/

	//Attempt #2
	//StressStraight segments[numNodes+1]; 
	for(int i = 0;i<DS.nodePositions.size()-1;i++)
	{

		int i2 = i+1;
		if(i2>DS.nodePositions.size()-1)
			i2 = 0;

		p0(0) = DS.nodePositions[i][0];
		p0(1) = DS.nodePositions[i][1];
		p0(2) = DS.nodePositions[i][2];

		//std::cout << "p0 " << p0 << std::endl;
		
		//Node 1 Positions
		p1(0) = DS.nodePositions[i2][0];
		p1(1) = DS.nodePositions[i2][1];
		p1(2) = DS.nodePositions[i2][2];

		//std::cout << "p1 " << p1 << std::endl;

		burgers(0) = DS.loops[0][0];
		burgers(1) = DS.loops[0][1];
		burgers(2) = DS.loops[0][2];

		//std::cout << "Calculating stress contribution from node " << i << " ("<< p0(0) <<"," << p0(1) << "," << p0(2) << ") to " << i2 << " ("<< p1(0) <<"," << p1(1) << "," << p1(2) << ")" << std::endl;

		tempMat = StressStraight(p0, p1, burgers).stress(positionVec);

		//std::cout << "Tempmat" << std::endl;
		//std::cout << tempMat << std::endl;

		//stressMat+=tempMat;

		stressMat=stressMat + tempMat;

		//std::cout << "Cumulative stressmat" << std::endl;
		//std::cout << stressMat << std::endl;
	}

	for(int i=0;i<3;i++)
		for(int j = 0; j<3; j++)
			appliedStress[i][j]=stressMat(i,j)*G;

	
	if( (useParametricStudy==1) && (useStress==1) ) //If performaing a parametric study with applied stress, factor in the applied stress.  
	{
		appliedStress[0][0] = appliedStress[0][0] + appliedPressure/3; 
		appliedStress[1][1] = appliedStress[1][1] + appliedPressure/3; 
		appliedStress[2][2] = appliedStress[2][2] + appliedPressure/3; 
	}

}

template <typename DislocationStructure>
void DSpositionToPeriodicStress(double pos[3], DislocationStructure DS,double appliedStress[][3])
{
	using namespace model;

	Eigen::Matrix<double,3,1> positionVec(pos[0]/b,pos[1]/b,pos[2]/b);

	Eigen::Matrix<double,3,3> stressMat(Eigen::Matrix<double,3,3>::Zero());

	Eigen::Matrix<double,3,3> tempMat;
	Eigen::Matrix<double,3,1> p0;  //Start of the dislocation segment
	Eigen::Matrix<double,3,1> p1;  //End of the dislocation segment
	Eigen::Matrix<double,3,1> burgers;


	for(int i = 0;i<DS.nodePositions.size()-1;i++)
	{

		int i2 = i+1;
		if(i2>DS.nodePositions.size()-1)
			i2 = 0;

		p0(0) = DS.nodePositions[i][0];
		p0(1) = DS.nodePositions[i][1];
		p0(2) = DS.nodePositions[i][2];

		//std::cout << "p0 " << p0 << std::endl;
		
		//Node 1 Positions
		p1(0) = DS.nodePositions[i2][0];
		p1(1) = DS.nodePositions[i2][1];
		p1(2) = DS.nodePositions[i2][2];

		//std::cout << "p1 " << p1 << std::endl;

		burgers(0) = DS.loops[0][0];
		burgers(1) = DS.loops[0][1];
		burgers(2) = DS.loops[0][2];

		//std::cout << "Calculating stress contribution from node " << i << " ("<< p0(0) <<"," << p0(1) << "," << p0(2) << ") to " << i2 << " ("<< p1(0) <<"," << p1(1) << "," << p1(2) << ")" << std::endl;

				

		tempMat = StressStraight(p0, p1, burgers).stress(positionVec);

		//std::cout << "Tempmat" << std::endl;
		//std::cout << tempMat << std::endl;


		stressMat=stressMat + tempMat; //Stress from original segment

		
		for(int j = 1;j< numImages;j++) //Now add a stress from the same segment, but repeated in images above and below the simulation box
		{
			p0(3) = p0(3) + j*L3;
			p1(3) = p1(3) + j*L3;

			stressMat = stressMat + StressStraight(p0, p1, burgers).stress(positionVec);

			p0(3) = p0(3) - j*L3;
			p1(3) = p1(3) - j*L3;

			stressMat = stressMat + StressStraight(p0, p1, burgers).stress(positionVec);
		}
		//std::cout << "Cumulative stressmat" << std::endl;
		//std::cout << stressMat << std::endl;
	}

	for(int i=0;i<3;i++)
		for(int j = 0; j<3; j++)
			appliedStress[i][j]=stressMat(i,j)*G;
}



void calcDirCosinesMat(double position[],double dirCosinesMat[][3])
//Function to calculate the direction cosines matric for a given stress tensor, which is calcualted from the vacancy position
{
	double appliedStress[3][3]; //Complex applied stress at original vac position
	double principleStress[3]; //Principle stresses at original vac position

	positionToAppliedStress(position, appliedStress); //Calculate applied stress given position

	princStress(appliedStress, principleStress);

	princDir(appliedStress, principleStress, dirCosinesMat);
}


void calcGradients(double position[], double moveGradient[], double delta)
//Functin to calculate the stress gradient, to define the vacancies movement direction
{

	//cout << "Begin gradient calculation" << endl;

	double appliedStress[3][3]; //Complex applied stress at original vac position
	double principleStress[3]; //Principle stresses at original vac position
	double dirCosinesMat[3][3]; //Principle directions of original point

	//Temp variables to exist during the loop
	double tempAppliedStress[3][3]; //Complex applied stress at original vac position
	double tempPrincipleStress[3]; //Array to hold principle stresses for gradient calc

	positionToAppliedStress(position, appliedStress); //Calculate applied stress given position

	princStress(appliedStress, principleStress);

	princDir(appliedStress, principleStress, dirCosinesMat);


	double deltaPos[3]; //Incremental change in either x,y,x directions

	for(int i=0; i<3;i++)
	{
		//Reset all the deltas to zero, except for the axis of interest
		deltaPos[0] = 0;
		deltaPos[1] = 0;
		deltaPos[2] = 0;

		deltaPos[i] = delta;


		//Convert to the axis basis of the principle directons basis
		//rotateToNewAxis(dirCosinesMat, deltaPos);
		rotateBackToOld(dirCosinesMat, deltaPos); //FLAGGED ADDED

		//cout << "Old Position: " << position[0] << " " << position[1] << " " << position[2] << endl;
		//Update positions and recalculate the principle stresses
		position[0] = position[0] + deltaPos[0];
		position[1] = position[1] + deltaPos[1];
		position[2] = position[2] + deltaPos[2];
		//cout << "New Position: " << position[0] << " " << position[1] << " " << position[2] << endl;


		positionToAppliedStress(position, tempAppliedStress); //Calculate applied stress given detla position
		princStress(tempAppliedStress, tempPrincipleStress);
		//print3x3(tempAppliedStress);
		//cout << endl;

		//Calculate the gradient
		moveGradient[i] = -((tempPrincipleStress[0]-principleStress[0]) 
						+ (tempPrincipleStress[1]-principleStress[1]) 
						+ (tempPrincipleStress[2]-principleStress[2]))/(3*delta);



		//Reset positions to the original position vector
		position[0] = position[0] - deltaPos[0];
		position[1] = position[1] - deltaPos[1];
		position[2] = position[2] - deltaPos[2];
	}


}

double calcOneGradient(double position[], double dirCosinesMat[][3], double deltaPos[])
//Functin to calculate the stress gradient, to define the vacancies movement direction
{

	double appliedStress[3][3]; //Complex applied stress at original vac position
	double principleStress[3]; //Principle stresses at original vac position
	//double dirCosinesMat[3][3]; //Principle directions of original point

	//Temp variables to exist during the loop
	double tempAppliedStress[3][3]; //Complex applied stress at original vac position
	double tempPrincipleStress[3]; //Array to hold principle stresses for gradient calc

	positionToAppliedStress(position, appliedStress); //Calculate applied stress given position

	princStress(appliedStress, principleStress);

	princDir(appliedStress, principleStress, dirCosinesMat);

	//Convert to the axis basis of the principle directons basis
	rotateToNewAxis(dirCosinesMat, deltaPos);

	//cout << "Old Position: " << position[0] << " " << position[1] << " " << position[2] << endl;
	//Update positions and recalculate the principle stresses
	position[0] = position[0] + deltaPos[0];
	position[1] = position[1] + deltaPos[1];
	position[2] = position[2] + deltaPos[2];
	//cout << "New Position: " << position[0] << " " << position[1] << " " << position[2] << endl;

	positionToAppliedStress(position, tempAppliedStress); //Calculate applied stress given detla position

	princStress(tempAppliedStress, tempPrincipleStress);
	//print3x3(tempAppliedStress);
	//cout << endl;


	//Calculate the gradient
	double driftGradient = ((tempPrincipleStress[0]-principleStress[0]) 
						+ (tempPrincipleStress[1]-principleStress[1]) 
						+ (tempPrincipleStress[2]-principleStress[2]))/(3*(deltaPos[0]+deltaPos[1]+deltaPos[2]));



	//Reset positions to the original position vector
	position[0] = position[0] - deltaPos[0];
	position[1] = position[1] - deltaPos[1];
	position[2] = position[2] - deltaPos[2];
	
	return driftGradient;
}


template <typename DislocationNetworkType>
void calc3x3x3Gradient(double position[], double allGradients[][3][3], double deltaPos, DislocationNetworkType& DN)
{
	double appliedStress[3][3]; //Complex applied stress at original vac position	

	//Temp variables to exist during the loop
	double tempAppliedStress[3][3]; //Complex applied stress at original vac position

	//positionToAppliedStress(position, appliedStress); //Calculate applied stress given position
	DDpositionToAppliedStress(position,DN, appliedStress);

	for(int deltaDirection=0;deltaDirection<3; deltaDirection++)
		{
		
		position[deltaDirection]= position[deltaDirection]+deltaPos;
		
		

		//std::cout << "Calculating stress for gradient..." << std::endl;
		DDpositionToAppliedStress(position,DN, tempAppliedStress); //Calculate the stress tensor at the new position
		
			
		
		for(int i = 0; i<3; i++)
			{
				for(int j = 0; j<3; j++)
				{
					allGradients[deltaDirection][i][j] = (appliedStress[i][j]-tempAppliedStress[i][j])/deltaPos; //Calculate the gradient
				}
			}
		
		position[deltaDirection]= position[deltaDirection]-deltaPos; //reset the position matrix

		}
}

template <typename DislocationStructure>
void calc3x3x3Gradient(double position[], double allGradients[][3][3], double deltaPos, DislocationStructure DS, int num)
{
	double appliedStress[3][3]; //Complex applied stress at original vac position	

	//Temp variables to exist during the loop
	double tempAppliedStress[3][3]; //Complex applied stress at original vac position

	//positionToAppliedStress(position, appliedStress); //Calculate applied stress given position
	DSpositionToAppliedStress(position,DS, appliedStress);

	for(int deltaDirection=0;deltaDirection<3; deltaDirection++)
		{
			
			position[deltaDirection]= position[deltaDirection]+deltaPos;
			
			//std::cout << "Calculating stress for gradient..." << std::endl;
			DSpositionToAppliedStress(position,DS, tempAppliedStress); //Calculate the stress tensor at the new position
			
			
			for(int i = 0; i<3; i++)
				{
					for(int j = 0; j<3; j++)
					{
						allGradients[deltaDirection][i][j] = (appliedStress[i][j]-tempAppliedStress[i][j])/deltaPos; //Calculate the gradient
					}
				}
			
			position[deltaDirection]= position[deltaDirection]-deltaPos; //reset the position matrix

		}
}

void calc3x3x3Gradient(double position[], double allGradients[][3][3], double deltaPos)
{
	double appliedStress[3][3]; //Complex applied stress at original vac position	

	//Temp variables to exist during the loop
	double tempAppliedStress[3][3]; //Complex applied stress at original vac position

	positionToAppliedStress(position, appliedStress); //Calculate applied stress given position
	//DDpositionToAppliedStress(position,DN, appliedStress);

	for(int deltaDirection=0;deltaDirection<3; deltaDirection++)
		{
		
		position[deltaDirection]= position[deltaDirection]+deltaPos;
				
		positionToAppliedStress(position, tempAppliedStress); //Calculate the stress tensor at the new position
		
		for(int i = 0; i<3; i++)
			{
				for(int j = 0; j<3; j++)
				{
					allGradients[deltaDirection][i][j] = (appliedStress[i][j]-tempAppliedStress[i][j])/deltaPos; //Calculate the gradient
				}
			}
		
		position[deltaDirection]= position[deltaDirection]-deltaPos; //reset the position matrix

		}
}

bool isInMesh(double vacPos[3])
{

	double posinB[3];

	for(int i=0;i<3;i++)
		posinB[i] = vacPos[i]/b;

	if (meshType==0)
		{
		if((posinB[0]<0) || (posinB[0]>L1) || (posinB[1]) < 0 || (posinB[1]>L2) || (posinB[2]<0) || (posinB[2]>L3))
			return false;
		else
			return true;
		}

	if (meshType==1)
		{
		if((pow(pow(posinB[0],2) + pow(posinB[1],2), 0.5)>L1) || (posinB[2]<0) || (posinB[2] >L2))
			return false;
		else
			return true;
		}
	
}

template <typename DislocationNetworkType>
void printstress(double pos[3], DislocationNetworkType& DN)
{

	using namespace model;

	//VectorDim positionVec(pos[0],pos[1],pos[2]);
	Eigen::Matrix<double,3,1> positionVec(pos[0],pos[1],pos[2]);
	//MatrixDim stressMat(MatrixDim::Zero());
	Eigen::Matrix<double,3,3> stressMat(Eigen::Matrix<double,3,3>::Zero());
	
	//stress(positionVec) = stressMat; 

	stressMat+=DN->stress(positionVec);

	std::cout << stressMat << std::endl;
}

void rewriteTemperature(double newTemp)
//Function to rewrite the W material input file with the updated temperature value
//Temperature is given in Kelvin - current function is designed to update the W material input file
{

	using namespace std;

	//Create a text file to hold the new W info called W_new, the delete the old W file and change the name of the new file

	ifstream InStream;
	InStream.open("../../MaterialsLibrary/W.txt"); //Original W file

	ofstream OutMaterialStream;
	OutMaterialStream.open("../../MaterialsLibrary/W_new.txt"); //W_new file to hold new W info

	string temp;
	string Indicator = "T="; //The last word that is written before the temperature line in the 	
	
	while (std::getline(InStream, temp) ) //Write the new base file to hold the W values
	{
		if (temp.find(Indicator) != std::string::npos)
		{

			OutMaterialStream << "T=" << newTemp << ";			# [K]		temperature";

			//InStream.ignore(256, '\n'); //Skip the next 2 /n markers
		}
		else
		{
			OutMaterialStream << temp;
		}

		OutMaterialStream << "\n";
	}

	OutMaterialStream.close();
	InStream.close(); 


	//Delete the old W material file
	remove("../../MaterialsLibrary/W.txt");

	//Rename the W_new file to have the name W
	char oldname[] = "../../MaterialsLibrary/W_new.txt";
	char newname[] = "../../MaterialsLibrary/W.txt";
	
	//Deletes the file if exists 
	if (rename(oldname, newname) != 0)
		perror("Error renaming material file\n");
	else
		std::cout << "Material file renamed successfully\n";
	
	
}

void resetRunningVariables()
//Function to call after running a parametric study trial to reset the running variables
{
	//////Reset the counters each parametric study trial!!//////
	RunningVacAbsorbed=0; 
	RunningLastVacNumber = 0; 
	RunningVacEmitted= 0;
	RunningLastVacEmitted=0; 
	RunningVacIDnum = 0; 
	totalGlobalTime = 0; 
	lastTotalGlobalTime = 0; 
	lastDistanceMoved = 0; 
	//////////////////////////////////////////

}

#endif
