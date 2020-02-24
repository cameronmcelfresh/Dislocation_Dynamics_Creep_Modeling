#include <iostream>
#include <cmath>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <string>
#include <iomanip>
#include <time.h>
#include <vector>
#include <complex>
#include <fstream>
#include <iomanip>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 
#include <DislocationStructure.h>
//#include <Vacancies_functions.h>


#ifndef vacancies_H_
#define vacancies_H_

using namespace std;
using namespace Eigen; 


class Vacancy
//Class to describe individual vacancies
{
public:
	Vacancy(); //empty constructor
	Vacancy(double positionValues[]); 
	Vacancy(double positionValues[], double thedriftGradients[]);

	void calcXjTimeStepGradient(); // Function to sample the jump distance from a gaussian, given a timestep - to be called when vacancies FIRST created

	template <typename DislocationNetworkType>
	void calcXjTimeStepGradient(DislocationNetworkType& DN); // Function to sample the jump distance from a gaussian, given a timestep

	//template <typename DislocationStructure>
	void calcXjTimeStepGradient(DislocationStructure& DS); // Function to sample the jump distance from a gaussian, given a timestep
	
	void calcJumpRates();

	void initializeVacancy(); // Function to call immediately after a vacancy is created - INITIALLY CALLED TO ALLOW VACANCIES TO BE GLOBAL

	template <typename DislocationNetworkType>
	void initializeVacancy(DislocationNetworkType& DN); // Function to call immediately after a vacancy is created OR immediately after it is moved

	//template <typename DislocationStructure>
	void initializeVacancy(DislocationStructure& DS); // Function to call immediately after a vacancy is created OR immediately after it is moved
	
	template <typename DislocationNetworkType>
	void moveVacancy(DislocationNetworkType& DN); // Function to call when a vacancy is selected to move - this function also updates the vacancies position

	//template <typename DislocationStructure>
	void moveVacancy(DislocationStructure& DS); // Function to call when a vacancy is selected to move - this function also updates the vacancies position

	//Variables inherent to the position/stress on vacancy
	int vacIDnum;
	double position[3];
	double driftGradients[3]; //Drift gradients are left in the PRINCIPLE ORIENTATIONS

	double dirCosinesMat[3][3]; //used to convert from principle orientations back to traditional axis

	//Arbitrary or calculated/intermediate variables
	double diffusionCoefficients[3]; //Dj
	double driftVelocity[3]; //ui, uk, uj

	//Statistically random variables
	double jumpDistances[3]; //Xij
	double jumpTimeStep[3]; //deltaTx, deltaTy, deltaTz
	double jumpRates[3]; //ri, rj, rk - rates of jumping in each principle direction
};

Vacancy::Vacancy()
//Empty constructor - drft distances and diffusion coefficients determined by random number generators
{
	//uniform_real_distribution<double> Uniformdistribution(0.0, 1.0);
	//random_device generator;
	//uniform_real_distribution<double> ZeroOnedistribution(0, 1); //Uniform distribution

	vacIDnum = RunningVacIDnum+1;
	RunningVacIDnum++;

	for(int i = 0; i<3; i++)
	{
		position[i] = (ZeroOnedistribution(generator)-0.5)*pow(10,-6);
		driftGradients[i] = ZeroOnedistribution(generator);
		diffusionCoefficients[i] = correlationFactor*nearestNeighbors*pow(b,2)*v*exp(-Eo/(Kb*Temp));
	}
	position[2]=0;	
	initializeVacancy();
}

Vacancy::Vacancy(double positionValues[])
{

	//uniform_real_distribution<double> Uniformdistribution(0.0, 1.0);
	//random_device generator;
	//uniform_real_distribution<double> ZeroOnedistribution(0, 1); //Uniform distribution

	vacIDnum = RunningVacIDnum+1;
	RunningVacIDnum++;

	for(int i = 0; i<3; i++)
	{
		position[i] = positionValues[i];
		driftGradients[i] = ZeroOnedistribution(generator);
		diffusionCoefficients[i] = correlationFactor*nearestNeighbors*pow(b,2)*v*exp(-Eo/(Kb*Temp));
	}
	//position[2]=0;	

	initializeVacancy();
}

Vacancy::Vacancy(double positionValues[], double thedriftGradients[])
{

	//uniform_real_distribution<double> Uniformdistribution(0.0, 1.0);

	vacIDnum = RunningVacIDnum+1;
	RunningVacIDnum++;

	for(int i = 0; i<3; i++)
	{
		driftGradients[i] = thedriftGradients[i];
		position[i] = positionValues[i];
		diffusionCoefficients[i] = correlationFactor*nearestNeighbors*pow(b,2)*v*exp(-Eo/(Kb*Temp));
	}
	position[2]=0;	

	initializeVacancy();
}

void Vacancy::calcXjTimeStepGradient()
//Function to sample (calculate) the Xj jump distance used in the rj calculation
{
	//random_device generator;
	//uniform_real_distribution<double> Stepdistribution((0.316 * pow(10,-10)), (0.316 * pow(10,-9))); //Uniform distribution for step distances

	for(int i=0; i<3; i++)
	{
		
		//This piece of code must remain until the vacancy is off of the z=0 plane
		if (i==2) //Set the jump distance and jump time for the z-step to the average of the 2, as a placeholder. 
			{
			jumpDistances[2] = (jumpDistances[0]+jumpDistances[1])/2;
			jumpTimeStep[2] = (jumpTimeStep[0]+jumpTimeStep[1])/2;
			continue;
			}
			
		

		double deltaT; //Guessed time step
		double deltaX; //Guessed step distance
		double xi;//Sampled step distance
		double u, STD; //Drift velocity and STD for sampling distribution
		double sampledTimeStep;

		double bestX, bestT;

		//Distribution to draw from for the initial deltaX step
		//uniform_real_distribution<double> deltaXdistribution(0, 1); 

		double error = 100;

		int iterations = 0;
		//for(int j= 0; j<10000; j++)
		//cout << "Searching for x[" << i << "] timestep..." << endl;
		do
		{

		//Guess step distance - do it iteratively to provide a randomly accepted value
		do{
		//deltaX = deltaXdistribution(generator); //Guessed X	
		deltaX = Stepdistribution(generator); //Guessed X from a uniform distribution akin to the unit cell length
		  } while(rand()%2>0);

		/*
		//Method 1
		double deltaPos[3] = {0,0,0};
		deltaPos[i] = deltaX;
		double driftGradient = calcOneGradient(position,  dirCosinesMat, deltaPos); //Calcualte the driftGradient for that given step size 
		cout << "driftGradient = " << driftGradient << endl;
		//Method 1 end
		*/

		/*
		//Method 2
		double tempMoveGradient[3];
		calcGradients(position,tempMoveGradient, deltaX);
		double driftGradient = tempMoveGradient[i];
		//cout << "Tempgradient = " << driftGradient << endl;
		//Method 2 end
		*/

		//METHOD 3 10/7/19 update
		double tempMoveGradientMatrix[3][3][3];
		
		calc3x3x3Gradient(position,tempMoveGradientMatrix, deltaX);

		double driftGradient = tempMoveGradientMatrix[i][0][0]+tempMoveGradientMatrix[i][1][1]+tempMoveGradientMatrix[i][2][2];
		//cout << "Tempgradient = " << driftGradient << endl;
		//METHOD 3 end
		

		//METHOD 1
		//Drift velocity caluclated from Marian notes
		//u = abs(2*deltaX*deltaX*v*Va/(Kb*Temp)*exp(-Eo/(Kb*Temp))*driftGradient); //Calculated Drift Velocity
		u = abs(2 * deltaX * deltaX * v * atomicVolume * (volumetricStrain/3) * exp(-Eo/(Kb*Temp)) * driftGradient / (Kb*Temp) ); //Calculated Drift Velocity

		//Guess time step
		double limit_num = log10(deltaX/u);	//Some fun math stuffs to try and find a more accurate range
		//METHOD 1 end

		//uniform_real_distribution<double> deltaTdistribution(pow(10, -10),pow(10, -13)); 

		uniform_real_distribution<double> deltaTdistribution(pow(10, limit_num-0.25),pow(10, limit_num+0.25)); 

		deltaT = deltaTdistribution(generator); //Guessed T

		STD = sqrt(2*diffusionCoefficients[i]*deltaX*deltaT); //Calcualte STD 

		normal_distribution<double> Testdistribution(u*deltaT, STD); 

		xi = Testdistribution(generator); //Sampled step distance

		
		//cout << "Iteration " << iterations << endl;
		//cout << "Diffusion Coefficient = " << diffusionCoefficients[0] << endl;
		//cout << "DriftGradient = " << driftGradient << endl;
		//cout << "u = " << u << endl;
		//cout << "Guessed X = " << deltaX << endl;
		//cout << "Guessed T = " << deltaT << endl;
		//cout << "Sampled Xi = " << xi << endl;
		//cout << (abs(xi-deltaX)/deltaX)*100<<  "%" << " error" << endl << endl;
		
		
		//If the error is less than error from previous trials, reassign the best values
		if ((abs(xi-deltaX)/deltaX)*100<error)
			{
			error = (abs(xi-deltaX)/deltaX)*100;
			bestX = deltaX;
			bestT = deltaT;
			}

		iterations++;

		if (error<kmcAcceptableError)
			driftGradients[i] = driftGradient; //If acceptable value, save the driftGradient Value 

		} while(error>kmcAcceptableError);

		//cout << "Best error @ " << error << "%" << endl;
		//cout << "Iterations: " << iterations << endl;


		//FLAGGED
		//If sucessful, update the jump distance and timestep

				
		if (driftGradients[i]<0)
			jumpDistances[i] = -deltaX;
		else
			jumpDistances[i] = deltaX;
		

		jumpTimeStep[i] = deltaT;

	}


	/*
	cout << "Best Results" << endl;
	cout << "Timesteps :" << jumpTimeStep[0] << "  " << jumpTimeStep[1] << "  " << jumpTimeStep[2] << endl;
	
	cout << "Jump Distances :" << jumpDistances[0] << "  " << jumpDistances[1] << "  " << jumpDistances[2] << endl;
	
	cout << "Position :" << position[0] << "  " << position[1] << "  " << position[2] << endl;
	cout << "driftGradients :" << driftGradients[0] << "  " << driftGradients[1] << "  " << driftGradients[2] << endl;
	*/

}

template <typename DislocationNetworkType>
void Vacancy::calcXjTimeStepGradient(DislocationNetworkType& DN)
//Function to sample (calculate) the Xj jump distance used in the rj calculation
{
	//random_device generator;
	//uniform_real_distribution<double> Stepdistribution((0.316 * pow(10,-10)), (0.316 * pow(10,-9))); //Uniform distribution for step distances

	for(int i=0; i<3; i++)
	{
		
		//This piece of code must remain until the vacancy is off of the z=0 plane
		//if (i==2) //Set the jump distance and jump time for the z-step to the average of the 2, as a placeholder. 
			//{
			//jumpDistances[2] = (jumpDistances[0]+jumpDistances[1])/2;
			//jumpTimeStep[2] = (jumpTimeStep[0]+jumpTimeStep[1])/2;
			//continue;
			//}
			
		

		double deltaT; //Guessed time step
		double deltaX; //Guessed step distance
		double xi;//Sampled step distance
		double u, STD; //Drift velocity and STD for sampling distribution
		double sampledTimeStep;

		double bestX, bestT;

		//Distribution to draw from for the initial deltaX step
		//uniform_real_distribution<double> deltaXdistribution(0, 1); 

		double error = 100;

		int iterations = 0;
		//for(int j= 0; j<10000; j++)
		//cout << "Searching for x[" << i << "] timestep..." << endl;
		do
		{

		//Guess step distance - do it iteratively to provide a randomly accepted value
		do{
		//deltaX = deltaXdistribution(generator); //Guessed X	
		deltaX = Stepdistribution(generator); //Guessed X from a uniform distribution akin to the unit cell length
		  } while(rand()%2>0);

		/*
		//Method 1
		double deltaPos[3] = {0,0,0};
		deltaPos[i] = deltaX;
		double driftGradient = calcOneGradient(position,  dirCosinesMat, deltaPos); //Calcualte the driftGradient for that given step size 
		cout << "driftGradient = " << driftGradient << endl;
		//Method 1 end
		*/

		/*
		//Method 2
		double tempMoveGradient[3];
		calcGradients(position,tempMoveGradient, deltaX);
		double driftGradient = tempMoveGradient[i];
		//cout << "Tempgradient = " << driftGradient << endl;
		//Method 2 end
		*/

		//METHOD 3 10/7/19 update
		double tempMoveGradientMatrix[3][3][3];
		
		//std::cout << "Calculating gradient" << std::endl;
		calc3x3x3Gradient(position,tempMoveGradientMatrix, deltaX, DN);

		double driftGradient = tempMoveGradientMatrix[i][0][0]+tempMoveGradientMatrix[i][1][1]+tempMoveGradientMatrix[i][2][2];
		//cout << "Tempgradient = " << driftGradient << endl;
		//METHOD 3 end
		

		//METHOD 1
		//Drift velocity caluclated from Marian notes
		//u = abs(2*deltaX*deltaX*v*Va/(Kb*Temp)*exp(-Eo/(Kb*Temp))*driftGradient); //Calculated Drift Velocity
		u = abs(2 * deltaX * deltaX * v * atomicVolume * (volumetricStrain/3) * exp(-Eo/(Kb*Temp)) * driftGradient / (Kb*Temp) ); //Calculated Drift Velocity

		//Guess time step
		double limit_num = log10(deltaX/u);	//Some fun math stuffs to try and find a more accurate range
		//METHOD 1 end

		//uniform_real_distribution<double> deltaTdistribution(pow(10, -10),pow(10, -13)); 

		uniform_real_distribution<double> deltaTdistribution(pow(10, limit_num-0.25),pow(10, limit_num+0.25)); 

		deltaT = deltaTdistribution(generator); //Guessed T

		STD = sqrt(2*diffusionCoefficients[i]*deltaX*deltaT); //Calcualte STD 

		normal_distribution<double> Testdistribution(u*deltaT, STD); 

		xi = Testdistribution(generator); //Sampled step distance

		
		//std::cout << "Iteration " << iterations << std::endl;
		//std::cout << "Diffusion Coefficient " << diffusionCoefficients[0] << endl;
		//std::cout << "DriftGradient = " << driftGradient << std::endl;
		//std::cout << "u = " << u << std::endl;
		//std::cout << "Guessed X = " << deltaX << std::endl;
		//std::cout << "Guseed T = " << deltaT << std::endl;
		//std::cout << "Sampled Xi = " << xi << std::endl;
		//std::cout << (abs(xi-deltaX)/deltaX)*100<<  "%" << " error" << std::endl << std::endl;
		
		//If the error is less than error from previous trials, reassign the best values
		if ((abs(xi-deltaX)/deltaX)*100<error)
			{
			error = (abs(xi-deltaX)/deltaX)*100;
			bestX = deltaX;
			bestT = deltaT;
			}

		iterations++;

		if (error<30)
			driftGradients[i] = driftGradient; //If acceptable value, save the driftGradient Value 

		} while(error>30);

		//cout << "Best error @ " << error << "%" << endl;
		//cout << "Iterations: " << iterations << endl;


		//FLAGGED
		//If sucessful, update the jump distance and timestep

				
		if (driftGradients[i]<0)
			jumpDistances[i] = -deltaX;
		else
			jumpDistances[i] = deltaX;
		

		jumpTimeStep[i] = deltaT;

	}


	/*
	cout << "Best Results" << endl;
	cout << "Timesteps :" << jumpTimeStep[0] << "  " << jumpTimeStep[1] << "  " << jumpTimeStep[2] << endl;
	
	cout << "Jump Distances :" << jumpDistances[0] << "  " << jumpDistances[1] << "  " << jumpDistances[2] << endl;
	
	cout << "Position :" << position[0] << "  " << position[1] << "  " << position[2] << endl;
	cout << "driftGradients :" << driftGradients[0] << "  " << driftGradients[1] << "  " << driftGradients[2] << endl;
	*/

}

//template <typename DislocationStructure>
void Vacancy::calcXjTimeStepGradient(DislocationStructure& DS)
//Function to sample (calculate) the Xj jump distance used in the rj calculation
{
	//random_device generator;
	//uniform_real_distribution<double> Stepdistribution((0.316 * pow(10,-10)), (0.316 * pow(10,-9))); //Uniform distribution for step distances

	for(int i=0; i<3; i++)
	{
		
		//This piece of code must remain until the vacancy is off of the z=0 plane
		//if (i==2) //Set the jump distance and jump time for the z-step to the average of the 2, as a placeholder. 
			//{
			//jumpDistances[2] = (jumpDistances[0]+jumpDistances[1])/2;
			//jumpTimeStep[2] = (jumpTimeStep[0]+jumpTimeStep[1])/2;
			//continue;
			//}
			
		

		double deltaT; //Guessed time step
		double deltaX; //Guessed step distance
		double xi;//Sampled step distance
		double u, STD; //Drift velocity and STD for sampling distribution
		double sampledTimeStep;

		double bestX, bestT;

		//Distribution to draw from for the initial deltaX step
		//uniform_real_distribution<double> deltaXdistribution(0, 1); 

		double error = 100;

		int iterations = 0;
		//for(int j= 0; j<10000; j++)
		//cout << "Searching for x[" << i << "] timestep..." << endl;
		do
		{

		//Guess step distance - do it iteratively to provide a randomly accepted value
		do{
		//deltaX = deltaXdistribution(generator); //Guessed X	
		deltaX = Stepdistribution(generator); //Guessed X from a uniform distribution akin to the unit cell length
		  } while(rand()%2>0);

		/*
		//Method 1
		double deltaPos[3] = {0,0,0};
		deltaPos[i] = deltaX;
		double driftGradient = calcOneGradient(position,  dirCosinesMat, deltaPos); //Calcualte the driftGradient for that given step size 
		cout << "driftGradient = " << driftGradient << endl;
		//Method 1 end
		*/

		/*
		//Method 2
		double tempMoveGradient[3];
		calcGradients(position,tempMoveGradient, deltaX);
		double driftGradient = tempMoveGradient[i];
		//cout << "Tempgradient = " << driftGradient << endl;

		//Method 2 end
		*/

		//METHOD 3 10/7/19 update
		double tempMoveGradientMatrix[3][3][3];
		
		//std::cout << "Calculating gradient" << std::endl;
		int placeholder = 1;
		calc3x3x3Gradient(position,tempMoveGradientMatrix, deltaX, DS,placeholder);

		double driftGradient = tempMoveGradientMatrix[i][0][0]+tempMoveGradientMatrix[i][1][1]+tempMoveGradientMatrix[i][2][2];
		//cout << "Tempgradient = " << driftGradient << endl;
		//METHOD 3 end
		

		//METHOD 1
		//Drift velocity caluclated from Marian notes
		//u = abs(2*deltaX*deltaX*v*Va/(Kb*Temp)*exp(-Eo/(Kb*Temp))*driftGradient); //Calculated Drift Velocity
		u = abs(2 * deltaX * deltaX * v * atomicVolume * (volumetricStrain/3) * exp(-Eo/(Kb*Temp)) * driftGradient / (Kb*Temp) ); //Calculated Drift Velocity

		//Guess time step
		double limit_num = log10(deltaX/u);	//Some fun math stuffs to try and find a more accurate range
		//METHOD 1 end

		//uniform_real_distribution<double> deltaTdistribution(pow(10, -10),pow(10, -13)); 

		uniform_real_distribution<double> deltaTdistribution(pow(10, limit_num-0.25),pow(10, limit_num+0.25)); 

		deltaT = deltaTdistribution(generator); //Guessed T

		STD = sqrt(2*diffusionCoefficients[i]*deltaX*deltaT); //Calcualte STD 


		normal_distribution<double> Testdistribution(u*deltaT, STD); 

		xi = Testdistribution(generator); //Sampled step distance

		
		//std::cout << "Iteration " << iterations << std::endl;
		//cout << "Diffusion Coefficient = " << diffusionCoefficients[0] << endl;
		//std::cout << "DriftGradient = " << driftGradient << std::endl;
		//std::cout << "u = " << u << std::endl;
		//std::cout << "Guessed X = " << deltaX << std::endl;
		//std::cout << "Guessed T = " << deltaT << std::endl;
		//std::cout << "Sampled Xi = " << xi << std::endl;
		//std::cout << (abs(xi-deltaX)/deltaX)*100<<  "%" << " error" << std::endl << std::endl;
		
		
		//If the error is less than error from previous trials, reassign the best values
		if ((abs(xi-deltaX)/deltaX)*100<error)
			{
			error = (abs(xi-deltaX)/deltaX)*100;
			bestX = deltaX;
			bestT = deltaT;
			}

		iterations++;

		if (error<30)
			driftGradients[i] = driftGradient; //If acceptable value, save the driftGradient Value 

		} while(error>30);

		//cout << "Best error @ " << error << "%" << endl;
		//cout << "Iterations: " << iterations << endl;


		//FLAGGED
		//If sucessful, update the jump distance and timestep

				
		if (driftGradients[i]<0)
			jumpDistances[i] = -deltaX;
		else
			jumpDistances[i] = deltaX;
		

		jumpTimeStep[i] = deltaT;

	}


	/*
	cout << "Best Results" << endl;
	cout << "Timesteps :" << jumpTimeStep[0] << "  " << jumpTimeStep[1] << "  " << jumpTimeStep[2] << endl;
	
	cout << "Jump Distances :" << jumpDistances[0] << "  " << jumpDistances[1] << "  " << jumpDistances[2] << endl;
	
	cout << "Position :" << position[0] << "  " << position[1] << "  " << position[2] << endl;
	cout << "driftGradients :" << driftGradients[0] << "  " << driftGradients[1] << "  " << driftGradients[2] << endl;
	*/

}

void Vacancy::calcJumpRates()
//Function to generate the ri , rj, rk, jump rates for individual vacancies
{
	for(int i = 0; i<3; i++)
		jumpRates[i] = diffusionCoefficients[i]/pow(jumpDistances[i],2);
}

void Vacancy::initializeVacancy()
//Function to reset the individual timesteps/jump distances when a vacancy is created or moved
{
	//std::cout << "Calculating dircosinesmat INITIAL..." << endl;
	calcDirCosinesMat(position, dirCosinesMat); //Recalcualte the direction cosines matrix

	/*
	for(int i = 0; i<3; i++)
		cout << "position[" << i << "] = " << position[i] << endl;
	*/
	//std::cout << "Calculating timestep INITIAL..." << endl;
	calcXjTimeStepGradient(); //Recalcualte jump distance, timesteps, and gradient values

	calcJumpRates(); //Recalcualte the individual timesteps
}

template <typename DislocationNetworkType>
void Vacancy::initializeVacancy(DislocationNetworkType& DN)
//Function to reset the individual timesteps/jump distances when a vacancy is created or moved
{
	//std::cout << "Calculating dircosinesmat considering dislocation network..." << endl;
	calcDirCosinesMat(position, dirCosinesMat); //Recalcualte the direction cosines matrix

	/*
	for(int i = 0; i<3; i++)
		cout << "position[" << i << "] = " << position[i] << endl;
	*/
	//std::cout << "Calculating timestep considering dislocation network..." << endl;
	calcXjTimeStepGradient(DN); //Recalcualte jump distance, timesteps, and gradient values

	calcJumpRates(); //Recalcualte the individual timesteps
}

//template <typename DislocationStructure>
void Vacancy::initializeVacancy(DislocationStructure& DS)
//Function to reset the individual timesteps/jump distances when a vacancy is created or moved
{
	//std::cout << "Calculating dircosinesmat considering dislocation network..." << endl;
	calcDirCosinesMat(position, dirCosinesMat); //Recalcualte the direction cosines matrix

	/*
	for(int i = 0; i<3; i++)
		cout << "position[" << i << "] = " << position[i] << endl;
	*/
	//std::cout << "Calculating timestep considering dislocation network..." << endl;
	calcXjTimeStepGradient(DS); //Recalcualte jump distance, timesteps, and gradient values

	calcJumpRates(); //Recalcualte the individual timesteps
}


template <typename DislocationNetworkType>
void Vacancy::moveVacancy(DislocationNetworkType& DN)
//Function to call when the vacancy is selected to move
{
	//First, augment the movegradient according to the appropriate deltaX,deltaY, deltaZ, and
	//then change the movegradient back to typical x,y,z orientations

	/*
	for(int i = 0; i<3;i++)
		driftGradients[i] = driftGradients[i]*jumpDistances[i];
	

	rotateBackToOld(dirCosinesMat, driftGradients); //rotate back to non-principle orientation
	*/

	double jumpIncrement[3];

	for(int i = 0;i<3; i++)
		{
		jumpIncrement[i] = jumpDistances[i];//[i]*jumpDistances[i]*driftGradients[i]; (2*jumpDistances[i]*jumpDistances[i]*v*Va*(1/(Kb*Temp))*exp(-Eo/(Kb*Temp))*driftGradients[i]);
		//cout << "JumpIncrement: " << jumpIncrement[i] << endl;
		}

	//rotateBackToOld(dirCosinesMat, jumpIncrement);
	//rotateToNewAxis(dirCosinesMat, jumpIncrement); //FLAGGED - changed because of how Eigen computes the direction cosines matrix

	//uppdate the positions given the calculated gradients
	position[0] = position[0] + jumpIncrement[0];
	position[1] = position[1] + jumpIncrement[1];
	position[2] = position[2] + jumpIncrement[2];

	/*
	cout << "Jump move distnace: " << endl;
	for(int i = 0;i<3; i++)
		cout << "[" << i << "] = " << jumpIncrement[i] << endl;
	*/
	initializeVacancy(DN); //Recalculating gradients, xij, rates, and timestep! 
}

//template <typename DislocationStructure>
void Vacancy::moveVacancy(DislocationStructure& DS)
//Function to call when the vacancy is selected to move
{
	//First, augment the movegradient according to the appropriate deltaX,deltaY, deltaZ, and
	//then change the movegradient back to typical x,y,z orientations

	/*
	for(int i = 0; i<3;i++)

		driftGradients[i] = driftGradients[i]*jumpDistances[i];
	

	rotateBackToOld(dirCosinesMat, driftGradients); //rotate back to non-principle orientation
	*/

	double jumpIncrement[3];

	for(int i = 0;i<3; i++)
		{
		jumpIncrement[i] = jumpDistances[i];//[i]*jumpDistances[i]*driftGradients[i]; (2*jumpDistances[i]*jumpDistances[i]*v*Va*(1/(Kb*Temp))*exp(-Eo/(Kb*Temp))*driftGradients[i]);
		//cout << "JumpIncrement: " << jumpIncrement[i] << endl;
		}

	//rotateBackToOld(dirCosinesMat, jumpIncrement);
	//rotateToNewAxis(dirCosinesMat, jumpIncrement); //FLAGGED - changed because of how Eigen computes the direction cosines matrix

	//uppdate the positions given the calculated gradients
	position[0] = position[0] + jumpIncrement[0];
	position[1] = position[1] + jumpIncrement[1];
	position[2] = position[2] + jumpIncrement[2];

	/*
	cout << "Jump move distnace: " << endl;
	for(int i = 0;i<3; i++)
		cout << "[" << i << "] = " << jumpIncrement[i] << endl;
	*/
	initializeVacancy(DS); //Recalculating gradients, xij, rates, and timestep! 
}


double calcRt(vector<Vacancy> &vacancyArray)
//Function to calculate the Rt given an array of vacancies 
{
	double Rt = 0;

	for(int i=0; i<vacancyArray.size(); i++)
	{
		for(int j = 0; j<3; j++)
		{
			Rt = Rt + vacancyArray[i].diffusionCoefficients[j]/(vacancyArray[i].jumpDistances[j]*vacancyArray[i].jumpDistances[j]);
		}
	}

	return Rt;
}

double calcGlobalT(vector<Vacancy> &vacancyArray)
//Function to calculate the delta timestep value given an Rt value
//where epsilon of -log(epsilon) is uniformly sampled from (0,1)
{
	//random_device generator;
	//uniform_real_distribution<double> ZeroOnedistribution(0, 1); //Uniform distribution

	double Rt = calcRt(vacancyArray); //Calculate the total Rt

	double deltaT;

	deltaT = -log(ZeroOnedistribution(generator))/Rt;

	return deltaT;
}

void runSimulation(vector<Vacancy> &vacancyArray, int steps)
//Function to run the vacancy simulation for a predeterminted number of steps
{
	//random_device generator;
	//uniform_real_distribution<double> ZeroOnedistribution(0, 1); //Uniform distribution

	int iters = 0;

	double globalTimeStep;
	double Rt;
	
	//Write out position data to a file
	ofstream OutXStream;
	OutXStream.open("vac_x_data.txt");

	ofstream OutYStream;
	OutYStream.open("vac_y_data.txt");

	ofstream OutPosStream;
	OutPosStream.open("vac_positions.txt");

	ofstream OutGlobalTimeStepStream;
	OutGlobalTimeStepStream.open("globalTimeStep.txt");

	//while(iters<steps)
	while(vacancyArray.size()>0)
	{

		Rt = calcRt(vacancyArray);

		double randomRate = Rt*ZeroOnedistribution(generator);

		int vacNum = -1;
		do
		{
			vacNum++;

			if ((vacancyArray[vacNum].position[0]==0) && (vacancyArray[vacNum].position[1]==0))
				continue;

			if (vacNum>(vacancyArray.size()-1))
				vacNum=0;

			for(int j =0;j<3;j++)
				{
					randomRate = randomRate-vacancyArray[vacNum].diffusionCoefficients[j]/pow(vacancyArray[vacNum].jumpDistances[j],2);
				}
		}while(randomRate>0);

		//cout << "Selected vacancy #" << vacNum << endl;

		if ((vacancyArray[vacNum].position[0]==0) && (vacancyArray[vacNum].position[1]==0)) //Failsafe
				continue;

		//vacancyArray[vacNum].moveVacancy(); 

		for(int i = 0; i<vacancyArray.size(); i++)
			{
			OutXStream << vacancyArray[i].position[0] << " ," ;
			OutYStream << vacancyArray[i].position[1] << " ," ;

				for(int j = 0;j<3; j++)
				{
					if(vacancyArray[i].position[j]<0)
						{
						OutPosStream << "  ";
						OutPosStream << std::scientific << setprecision(15) << std::right << vacancyArray[i].position[j];
						}
					else
						{
						OutPosStream << "   ";
						OutPosStream << std::scientific << setprecision(15) << std::right << vacancyArray[i].position[j];
						}
				}
			}

		OutXStream << endl;
		OutYStream << endl;
		OutPosStream << endl;


		//Check to see if the vacancy has "reached" the dislcoation - if so, remove it and replace it with a new
		//randomly placed vacancy
		if (sqrt(pow(vacancyArray[vacNum].position[0],2)+pow(vacancyArray[vacNum].position[1],2))<5*pow(10,-8))
			{
			//vacancyArray.erase(vacancyArray.begin()+vacNum);
			vacancyArray[vacNum].position[0] = 0;
			vacancyArray[vacNum].position[1] = 0;

/*
			if (vacancyArray.size()==0) //Exit out of the simulation if all the vacancies have reached the dislocation
				{
				cout << "Done!" << endl;
				return;
				}
*/

			int exitVal = 1;
			for (int i = 0; i<vacancyArray.size();i++) //Exit out of the simulation if all the vacancies have reached the dislocation (0,0)
				{
				if((vacancyArray[i].position[0]!=0) || (vacancyArray[i].position[0]!=0))
					{
						exitVal = 0;
						break;
					}					
				}

				if(exitVal==1)
				{
					globalTimeStep = calcGlobalT(vacancyArray); 
					OutGlobalTimeStepStream << globalTimeStep << endl;

					OutXStream.close();
					OutYStream.close();
					OutPosStream.close();
					OutGlobalTimeStepStream.close();
					std::cout << "All vacancies gone!" << endl;
					return;
				}
			}


		globalTimeStep = calcGlobalT(vacancyArray); //Calculate the global timestep
		OutGlobalTimeStepStream << globalTimeStep << endl;


		iters++;
	}
	
	OutXStream.close();
	OutYStream.close();
	OutPosStream.close();
	OutGlobalTimeStepStream.close();
}

void vacRecount(vector<Vacancy> &vacancyArray, DislocationStructure DisStruct, double timeStep)
//Function to recount the number of vacancies and add or remove vacancies to keep the original vacancy concentration
//This function assists in nullifying any effects that may occur from change in the vacancy number due to absorption or emission
{

	std::cout <<"Checking total vacancies to ensure desired #..." << std::endl;

	////Calculate the total absorption and emission rates to find the balancing rate////

	double absorptionRate = RunningVacAbsorbed/totalGlobalTime; //Total absorption rate
	double emissionRate = RunningVacEmitted/totalGlobalTime; //Total emission rate

	double addingRate = (absorptionRate-emissionRate); //The rate of vacancy adding (+) or subtracting (-) that must occur to keep a constant concentration

	double prob = abs(addingRate*timeStep); //Probability of adding/removing a vacancy given the current timestep

	double num_of_vacs = static_cast<double>(vacancyArray.size()); //Length of the vacancy vector aka # of vacancies

	std::cout << "The balancing rate for const. concentration is --> " << addingRate << " vacs/sec --> " << prob*100 + 100*(vacNum-num_of_vacs)*10/vacNum << "%" << " [adjusted]" << std::endl;

	std::normal_distribution<double> VacsToAdd(addingRate*timeStep + (vacNum-num_of_vacs)*10/vacNum, pow(abs(addingRate*timeStep), 0.5)); //Noraml distribution to determine how many vacancies to add/subtract

	double incrementalVacNum = VacsToAdd(generator); //Number of vacancies to add or subtract during this iteration

	std::cout << "Randomly generated incrementalVacNum = " << incrementalVacNum << std::endl;

	incrementalVacNum = round( incrementalVacNum ); //Round to the closest integer number after adding a percentage change

	if(abs(incrementalVacNum)>0)
		{
	
		std::cout << "Sampled additional vacancy number: " << incrementalVacNum << std::endl;

		while(incrementalVacNum>0) //If there are too few vacancies, add more until the constantVacNum is reached
			{
			std::cout << vacancyArray.size() << " vacancies, adding vacancy" << std::endl;

			////////Add a vacancy somewhere in the mesh and initialize it////////
			std::uniform_real_distribution<double> L1Dist(-L1/2,L1/2);
			std::uniform_real_distribution<double> L2Dist(-L3/2,L2/2);
			std::uniform_real_distribution<double> L3Dist(-L3/2,L3/2);

			double pos[3]; 

			pos[0] = L1Dist(generator);	
			pos[1] = L2Dist(generator);	
			pos[2] = L3Dist(generator);	

			Vacancy newVacancy;//Create	

			vacancyArray.insert(vacancyArray.begin(),newVacancy ); //Add a vacancy to the current vector of vacancies

			for(int i = 0; i<3; i++)
				vacancyArray[0].position[i] = pos[i]*b;

			vacancyArray[0].initializeVacancy(DisStruct); //initialize the vacancy!
			
			RunningVacIDnum++;
			
			vacancyArray[0].vacIDnum=RunningVacIDnum;

			incrementalVacNum--; 


			RunningBalancedVacs++; //Increase the running counter
			}


		while(incrementalVacNum<0 && vacancyArray.size()>1) //If there are too many vacancies, remove until the constantVacNum is reached
			{
			std::cout << vacancyArray.size() << " vacancies, removing vacancy" << std::endl;

			std::uniform_int_distribution<> vacToRemove(0,vacancyArray.size()-1); //Uniform distribution to select a random vacancy to remove

			int removedVac = vacToRemove(generator); //random vacancy number to remove

			vacancyArray.erase(vacancyArray.begin()+removedVac); //Remove the vacancy 

			incrementalVacNum++; 

			RunningBalancedVacs--; //Decrease the running counter
			}

		}


	std::cout <<"Total vacancies added/subtracted = " << RunningBalancedVacs << std::endl;

}

template <typename DislocationNetworkType>
void singleVacancyEvent(vector<Vacancy> &vacancyArray, DislocationNetworkType& DN, long int runIDnum, DislocationStructure DisStruct)
//Function to run the vacancy simulation step-wise
{
	//random_device generator;
	//uniform_real_distribution<double> ZeroOnedistribution(0, 1); //Uniform distribution


	double globalTimeStep;
	double Rt;
	

	if(vacancyArray.size()>0)
	{

		Rt = calcRt(vacancyArray);

		double randomRate = Rt*ZeroOnedistribution(generator); //Select a random vacancy to move

		int vacNum = -1;
		do
		{
			vacNum++;

			if ((vacancyArray[vacNum].position[0]==0) && (vacancyArray[vacNum].position[1]==0))
				continue;

			if (vacNum>(vacancyArray.size()-1))
				vacNum=0;

			for(int j =0;j<3;j++)
				{
					randomRate = randomRate-vacancyArray[vacNum].diffusionCoefficients[j]/pow(vacancyArray[vacNum].jumpDistances[j],2);
				}
		}while(randomRate>0); //Move through the vacancy list to find the proper timsetp

		std::cout << "Selected vacancy #" << vacancyArray[vacNum].vacIDnum << " to move" << std::endl;


		vacancyArray[vacNum].moveVacancy(DN); //Move the vacancy!

		//////Check to see if the vacancy has exited the crystal - if so, reset its position on the other side to maintain perioidc boundary conditions////////
		int moved = 0;

			if (vacancyArray[vacNum].position[0]>L1*b/2)
				{				
				vacancyArray[vacNum].position[0] = vacancyArray[vacNum].position[0] - L1*b;
				moved++;
				}
			if (vacancyArray[vacNum].position[0]<-L1*b/2)
				{				
				vacancyArray[vacNum].position[0] = vacancyArray[vacNum].position[0] + L1*b;
				moved++;
				}

			if (vacancyArray[vacNum].position[1]>L2*b/2)
				{				
				vacancyArray[vacNum].position[1] = vacancyArray[vacNum].position[1] - L2*b;
				moved++;
				}

			if (vacancyArray[vacNum].position[1]<-L2*b/2)
				{				
				vacancyArray[vacNum].position[1] = vacancyArray[vacNum].position[1] + L2*b;
				moved++;
				}

			if (vacancyArray[vacNum].position[2]>L3*b/2)
				{				
				vacancyArray[vacNum].position[2] = vacancyArray[vacNum].position[2] - L3*b;
				moved++;
				}

			if (vacancyArray[vacNum].position[2]<-L3*b/2)
				{				
				vacancyArray[vacNum].position[2] = vacancyArray[vacNum].position[2] + L3*b;	
				moved++;
				}

		if(moved>0)			
			vacancyArray[vacNum].initializeVacancy(DN); //Reinitialize the vacancy if it has been placed into a new position in the simulation box
		//////////////////////////////////////////////////

		//Write out position data to a file if outputV ==1
		if(outputV==1)
		{
			ofstream OutXStream;
			OutXStream.open("./V/vac_x_data.txt", std::ios_base::app);

			ofstream OutYStream;
			OutYStream.open("./V/vac_y_data.txt", std::ios_base::app);

			std::string fileString= "./V/V_" + std::to_string(runIDnum) + ".txt";  //Creates file name based on the current runID
			ofstream OutPosStream;
			OutPosStream.open(fileString);

			for(int i = 0; i<vacancyArray.size(); i++)
				{
				OutXStream << vacancyArray[i].position[0] << " ," ;
				OutYStream << vacancyArray[i].position[1] << " ," ;


				OutPosStream << vacancyArray[i].vacIDnum << "  ";

					for(int j = 0;j<3; j++)
					{
						
						if(vacancyArray[i].position[j]<0)
							{
							OutPosStream << "  ";
							OutPosStream << std::scientific << setprecision(15) << std::right << vacancyArray[i].position[j]/b;
							}
						else
							{
							OutPosStream << "   ";
							OutPosStream << std::scientific << setprecision(15) << std::right << vacancyArray[i].position[j]/b;
							}

					}

					OutPosStream << endl;

				}

				OutXStream << endl;
				OutYStream << endl;

				OutXStream.close();
				OutYStream.close();
				OutPosStream.close();

		}
		


		//Calculate the and print the global timestep
		globalTimeStep = calcGlobalT(vacancyArray); 

		totalGlobalTime=totalGlobalTime + globalTimeStep;

		if(outputGlobalTimeStep==1) //Output the global timestep 
		{
			ofstream OutGlobalTimeStepStream;
			OutGlobalTimeStepStream.open("./V/globalTimeStep.txt",std::ios_base::app);

			OutGlobalTimeStepStream << globalTimeStep << endl;
			OutGlobalTimeStepStream.close();
		}
	}

	
}


void singleVacancyEvent(vector<Vacancy> &vacancyArray, long int runIDnum, DislocationStructure& DisStruct)
//Function to run the vacancy simulation step-wise
{
	//random_device generator;
	//uniform_real_distribution<double> ZeroOnedistribution(0, 1); //Uniform distribution


	double globalTimeStep;
	double Rt;
	

	if(vacancyArray.size()>0)
	{

		Rt = calcRt(vacancyArray);

		double randomRate = Rt*ZeroOnedistribution(generator);

		int vacNum = -1;
		do
		{
			vacNum++;

			if ((vacancyArray[vacNum].position[0]==0) && (vacancyArray[vacNum].position[1]==0))
				continue;

			if (vacNum>(vacancyArray.size()-1))
				vacNum=0;

			for(int j =0;j<3;j++)
				{
					randomRate = randomRate-vacancyArray[vacNum].diffusionCoefficients[j]/pow(vacancyArray[vacNum].jumpDistances[j],2);
				}
		}while(randomRate>0);

		std::cout << "Selected vacancy #" << vacancyArray[vacNum].vacIDnum << " to move" << std::endl;


		vacancyArray[vacNum].moveVacancy(DisStruct); //Move the vacancy!

		//Check to see if the vacancy has exited the crystal - if so, reset its position on the other side to maintain perioidc boundary conditions
		int moved = 0;

			if (vacancyArray[vacNum].position[0]>L1*b/2)
				{				
				vacancyArray[vacNum].position[0] = vacancyArray[vacNum].position[0] - L1*b;
				moved++;
				}
			if (vacancyArray[vacNum].position[0]<-L1*b/2)
				{				
				vacancyArray[vacNum].position[0] = vacancyArray[vacNum].position[0] + L1*b;
				moved++;
				}

			if (vacancyArray[vacNum].position[1]>L2*b/2)
				{				
				vacancyArray[vacNum].position[1] = vacancyArray[vacNum].position[1] - L2*b;
				moved++;
				}

			if (vacancyArray[vacNum].position[1]<-L2*b/2)
				{				
				vacancyArray[vacNum].position[1] = vacancyArray[vacNum].position[1] + L2*b;
				moved++;
				}

			if (vacancyArray[vacNum].position[2]>L3*b/2)
				{				
				vacancyArray[vacNum].position[2] = vacancyArray[vacNum].position[2] - L3*b;
				moved++;
				}

			if (vacancyArray[vacNum].position[2]<-L3*b/2)
				{				
				vacancyArray[vacNum].position[2] = vacancyArray[vacNum].position[2] + L3*b;	
				moved++;
				}

		if(moved>0)			
			vacancyArray[vacNum].initializeVacancy(DisStruct);


		//Write out position data to a file IF outputV ==1

		if(outputV==1)
		{
			ofstream OutXStream;
			OutXStream.open("./V/vac_x_data.txt", std::ios_base::app);

			ofstream OutYStream;
			OutYStream.open("./V/vac_y_data.txt", std::ios_base::app);

			std::string fileString= "./V/V_" + std::to_string(runIDnum) + ".txt";  //Creates file name based on the current runID
			ofstream OutPosStream;
			OutPosStream.open(fileString);

			for(int i = 0; i<vacancyArray.size(); i++)
				{
				OutXStream << vacancyArray[i].position[0] << " ," ;
				OutYStream << vacancyArray[i].position[1] << " ," ;


				OutPosStream << vacancyArray[i].vacIDnum << "  ";

					for(int j = 0;j<3; j++)
					{
						
						if(vacancyArray[i].position[j]<0)
							{
							OutPosStream << "  ";
							OutPosStream << std::scientific << setprecision(15) << std::right << vacancyArray[i].position[j]/b;
							}
						else
							{
							OutPosStream << "   ";
							OutPosStream << std::scientific << setprecision(15) << std::right << vacancyArray[i].position[j]/b;
							}

					}

					OutPosStream << endl;

				}

				OutXStream << endl;
				OutYStream << endl;

				OutXStream.close();
				OutYStream.close();
				OutPosStream.close();

		}
		


		//////Calculate the and print the global timestep//////
		globalTimeStep = calcGlobalT(vacancyArray); 

		totalGlobalTime=totalGlobalTime + globalTimeStep;

		if(outputGlobalTimeStep==1) //Output the global timestep 
		{
			ofstream OutGlobalTimeStepStream;
			OutGlobalTimeStepStream.open("./V/globalTimeStep.txt",std::ios_base::app);

			OutGlobalTimeStepStream << globalTimeStep << endl;
			OutGlobalTimeStepStream.close();
		}
	}


	/////Go through vacancy emission procedure if useEmission is turned on////
	if(useEmission==1)
		vacancyEmission(vacancyArray, DisStruct, globalTimeStep);

	/////Add or subtract vacancies to keep a constant vacancy concentration if preset////
	if(requireConstantVacs==1)
		vacRecount(vacancyArray,  DisStruct, globalTimeStep);
		

}


template <typename DislocationNetworkType>
void findVacancyIntersections(vector<Vacancy> &vacancyArray, DislocationStructure& DisStruct, DislocationNetworkType& DN)
//Function to find all the absorptions of the vacancies to segments, then update the nodal positons
{
	int vacsAbsorbed=0;

	double distanceInfo[3];

	for(int i = 0; i<vacancyArray.size(); i++)
			{
				findClosestSegment(vacancyArray[i].position, DisStruct,  distanceInfo);
				// distanceInfo[0] = distance to segment in b units
				// distanceInfo[1] = nodeID1
				// distanceInfo[2] = nodeID2
				//std::cout << "Vac " << vacancyArray[i].vacIDnum << " closest segment is " << distanceInfo[0] << " b" << std::endl; 

				//Skip the segment if it is a screw dislocation aka a connecting segment
				//Confirm by seeing if they are on the same z-position
				if( abs( DisStruct.nodePositions[distanceInfo[1]][2] - DisStruct.nodePositions[distanceInfo[2]][2] ) ==0 ) 
					continue;


				if (distanceInfo[0]<distToAbsorbption)
				//if (pow(pow(vacancyArray[i].position[0]/b, 2) +pow(vacancyArray[i].position[1]/b, 2),0.5) < distToAbsorbption)
					{

					/////////////////////////////
					//Calculate the climb height
					/////////////////////////////
					double atomicVolumeinb = atomicVolume*pow(b,-3); //Atomic volume in b units

					double segmentLength = pow(pow(DisStruct.nodePositions[distanceInfo[1]][0]-DisStruct.nodePositions[distanceInfo[2]][0],2) + 
									pow(DisStruct.nodePositions[distanceInfo[1]][1]-DisStruct.nodePositions[distanceInfo[2]][1],2) +
									pow(DisStruct.nodePositions[distanceInfo[1]][2]-DisStruct.nodePositions[distanceInfo[2]][2],2), 0.5);

					double h = ( atomicVolumeinb*(1+volumetricStrain)/ segmentLength); //Climb height in b units according to volume swept out by absorption of one vacancy

					//Get the plane normal -- begin by getting the loop number
					int loopIDnum;

					for( auto const& [key, val] : DisStruct.nodeConnectivity)
							{
							if(key[0]==distanceInfo[0] || key[1]==distanceInfo[0] || key[0]==distanceInfo[1] || key[1]==distanceInfo[1])
								loopIDnum= val[0];
							}

					double planeNormal[3]; //Plane normal in unit vector form

					for( auto const& [key, val] : DisStruct.loops)
							{
							if(key==loopIDnum)
								{
									planeNormal[0] = val[3];
									planeNormal[1] = val[4];
									planeNormal[2] = val[5]; //Get the plane normals from the correct loop ID

									break;
								}
							}

					/////////////////////////////////////////////////////////////////////////////////////////
					//Execute vacancy absorption and place the vacancy in a new position	
					/////////////////////////////////////////////////////////////////////////////////////////				
					
					vacsAbsorbed++;
					RunningVacAbsorbed++;
					
					if(useEmission==0)
					{
						//Absorb the vacancy and randomly place it somewhere else
						double pos[3];

						std::uniform_real_distribution<double> L1Dist(-L1/2,L1/2);
						std::uniform_real_distribution<double> L2Dist(-L3/2,L2/2);
						std::uniform_real_distribution<double> L3Dist(-L3/2,L3/2);

						pos[0] = L1Dist(generator);	
						pos[1] = L2Dist(generator);	
						pos[2] = L3Dist(generator);	

						for(int j = 0; j<3;j++)
							vacancyArray[i].position[j] = pos[j]*b;

						vacancyArray[i].initializeVacancy(DN);

						RunningVacIDnum++;
						vacancyArray[i].vacIDnum=RunningVacIDnum;
					}
					else
					{
						vacancyArray.erase(vacancyArray.begin()+i); //Absorb and remove the vacancy if vacancy emission is turned on
					}

					/////////////////////////////////////////////////////////////////////////////////////////				
					//Execute CLIMB of the segment that the vacancy was absorbed into - via movement of the two segments
					/////////////////////////////////////////////////////////////////////////////////////////				

					Eigen::Matrix<double,3,1> newNodePosition(0,0,0);			

					    for (auto& nodeIter : DN->nodes())
					    {
						if (nodeIter.first==distanceInfo[1])
							{
							DisStruct.nodePositions[distanceInfo[1]][0]+=h*planeNormal[0];
							DisStruct.nodePositions[distanceInfo[1]][1]+=h*planeNormal[1];
							DisStruct.nodePositions[distanceInfo[1]][2]+=h*planeNormal[2];

							newNodePosition(0) = DisStruct.nodePositions[distanceInfo[1]][0];
							newNodePosition(1) = DisStruct.nodePositions[distanceInfo[1]][1];
							newNodePosition(2) = DisStruct.nodePositions[distanceInfo[1]][2];

							nodeIter.second->set_P(newNodePosition);
							}

						if (nodeIter.first==distanceInfo[2])
							{

							DisStruct.nodePositions[distanceInfo[2]][0]+=h*planeNormal[0];
							DisStruct.nodePositions[distanceInfo[2]][1]+=h*planeNormal[1];
							DisStruct.nodePositions[distanceInfo[2]][2]+=h*planeNormal[2];

							newNodePosition(0) = DisStruct.nodePositions[distanceInfo[2]][0];
							newNodePosition(1) = DisStruct.nodePositions[distanceInfo[2]][1];
							newNodePosition(2) = DisStruct.nodePositions[distanceInfo[2]][2];

							nodeIter.second->set_P(newNodePosition);
							}
					    } 

					}
			}

		std::cout << vacsAbsorbed << " vacancies absorbed" << std::endl;
		std::cout << RunningVacAbsorbed << " total vacancies absorbed" << std::endl;
}

void remesh(DislocationStructure& DisStruct)
//Function to cycle through the mesh and add or remove any nodes as necessary
{
	int isMeshed = 0;


	while(isMeshed!=1)
	{



		isMeshed = 1; //Assume completely meshed, then go through each node-node segment

		for( auto const& [key, val] : DisStruct.nodeConnectivity)
		{

			if( (key[0] ==numNodes+1) || (key[0] ==numNodes) || (key[0] ==numNodes-1) )
			{
			continue; //To skip all the boundary segments which shouldn't be remeshed
			}

			//std::cout << "Node " << key[0] << " and " << key[1] << std::endl;

			double x_dif = DisStruct.nodePositions[key[0]][0] - DisStruct.nodePositions[key[1]][0];  
			double y_dif = DisStruct.nodePositions[key[0]][1] - DisStruct.nodePositions[key[1]][1];
			double z_dif = DisStruct.nodePositions[key[0]][2] - DisStruct.nodePositions[key[1]][2];

			x_dif = abs(x_dif);
			y_dif = abs(y_dif);
			z_dif = abs(z_dif);


			int firstNode = key[0];
			int secondNode = key[1];

			//std::cout << "x_dif = " << x_dif << std::endl;
			//std::cout << "y_dif = " << y_dif << std::endl;
			//std::cout << "z_dif = " << z_dif << std::endl;


			//First check to see if there is an angled segment, if so, add a node to make it straight
			if (x_dif>0 && z_dif>0)
			{

				if(key[0]==key[1])
					{
					continue; //Sanity Check
					}

				//std::cout << "Looping at curvature node" << std::endl;

				int newNodeNum=0;

				//Find the biggest nodeID number and add one to it
				for( auto const& [nodeKeys, nodeVals] : DisStruct.nodeConnectivity)
				{
					if(nodeKeys[1]>newNodeNum) 
						{				
							newNodeNum=nodeKeys[1];
						}	

					if(nodeKeys[0]>newNodeNum) 
						{				
							newNodeNum=nodeKeys[0];
						}
				}
	

				newNodeNum = newNodeNum+1; //Make your node num value one higher than the largest value

				//std::cout << "(1)New node num = " << newNodeNum << std::endl;
				//std::cout << "(1)key[1] = " << secondNode << std::endl;

				//Create and insert the new node
				std::array<double,6> position;

				
				//Select the lowest x-value of the two
				position[0]= min(DisStruct.nodePositions[key[0]][0], DisStruct.nodePositions[key[1]][0]);
				position[1]= 0;

				//Select the furthest away z node to provide the x position
				if (position[0] ==DisStruct.nodePositions[key[0]][0])
					{
					position[2] = DisStruct.nodePositions[key[1]][2];
					}			
				else
					{
					position[2] = DisStruct.nodePositions[key[0]][2];
					}

				position[3] = 0;
				position[4] = 0;
				position[5] = 0;

				//std::cout << "Curvature between node "<< key[0] << " at :" << DisStruct.nodePositions[key[0]][0] << ", " << DisStruct.nodePositions[key[0]][1] << ", " << DisStruct.nodePositions[key[0]][2] << std::endl;
				//std::cout << " and node " << key[1] << " at " << DisStruct.nodePositions[key[1]][0] << ", " << DisStruct.nodePositions[key[1]][1] << ", " << DisStruct.nodePositions[key[1]][2] << std::endl;
				//std::cout << "Adding node " << newNodeNum << " at : " << position[0] << ", " << position[1] << ", " << position[2] << std::endl;

				//std::cout << "(2)New node num = " << newNodeNum << std::endl;
				//std::cout << "(2)key[1] = " << secondNode << std::endl;

				//std::cout << "Inserting node " << newNodeNum << std::endl;
				DisStruct.nodePositions.insert( std::make_pair(newNodeNum, position)  );

				//Update the node connectivity 

				//Delete old connectivity
				std::array<int, 2> nodes;

				nodes[0] = key[0];
				nodes[1] = key[1];

				//std::cout << "Deleting segment " << nodes[0] <<"->" <<nodes[1] <<std::endl;
				DisStruct.nodeConnectivity.erase(nodes); //erase the old connectivity value

				//Insert the new connectivity values
				std::array<int,2> connectivity;

				connectivity[0] = firstNode;
				connectivity[1] = newNodeNum;

				//std::cout << "(3)New node num = " << newNodeNum << std::endl;
				//std::cout << "(3)key[1] = " << secondNode << std::endl;

				assert(key[0]!=newNodeNum); //Sanity Check

				std::array<int,2> nodeInfo;
				nodeInfo[0] = 0;
				nodeInfo[1] = 0;

				//std::cout << "Adding segment " << connectivity[0] <<"->" <<connectivity[1] <<std::endl;
				DisStruct.nodeConnectivity.insert( std::make_pair(connectivity, nodeInfo)  );

				//Second connectivity value
				std::array<int,2> connectivity1;
				connectivity1[0] = newNodeNum;
				connectivity1[1] = secondNode;

				//std::cout << "(4)New node num = " << newNodeNum << std::endl;
				//std::cout << "(4)key[1] = " << secondNode << std::endl;

				assert(secondNode!=newNodeNum); //Sanity Check

				//std::cout << "Adding segment " << connectivity1[0] <<"->" <<connectivity1[1] <<std::endl;
				DisStruct.nodeConnectivity.insert( std::make_pair(connectivity1, nodeInfo)  );


				assert(DisStruct.nodePositions.size()==DisStruct.nodeConnectivity.size()); //Make sure we maintain the required lengths

				isMeshed = 0; //Do not allow the remeshing to end
				break; //Break out and repeat the loop if remeshing was needed

			}


			//Remove a node if two nodes now occupy the same position
			if((x_dif==0) && (y_dif==0) && (z_dif==0) && (isMeshed==1))
			{
				
				if(key[0]==key[1])
					{
					continue; //Sanity Check
					}

				//std::cout << "Overlapping nodes found" << std::endl;

				//std::cout << "Node " << key[0] << " :";
				//for(int i = 0; i <3; i++)
				//	{
				//	std::cout << DisStruct.nodePositions[key[0]][i];
				//	}
				
				//std::cout << std::endl;

				//std::cout << "Node " << key[1] << " :";
				//for(int i = 0; i <3; i++)
				//	{
				//	std::cout << DisStruct.nodePositions[key[1]][i];
				//	}
										
				//std::cout << std::endl;

				int nextNode;			

				//First, figure out the node connecting to the second node
				for( auto const& [key1, val1] : DisStruct.nodeConnectivity)
				{
					if(key1[0]==key[1]) 
						{				
							nextNode = key1[1]; //Finding the next node in line
						}	

				}

				//Remove the old node
				DisStruct.nodePositions.erase(key[1]);

				//Update the node connectivity//

				//Delete old connectivity
				std::array<int, 2> nodes;

				nodes[0] = firstNode;
				nodes[1] = secondNode;

				DisStruct.nodeConnectivity.erase(nodes); //erase the old connectivity value

				nodes[0] = secondNode;
				nodes[1] = nextNode;

				DisStruct.nodeConnectivity.erase(nodes); //erase the old connectivity value

				//Insert the new connectivity values
				std::array<int,2> connectivity;
				connectivity[0] = firstNode;
				connectivity[1] = nextNode;

				assert(key[0]!=nextNode); //Sanity Check

				std::array<int,2> nodeInfo;
				nodeInfo[0] = 0;
				nodeInfo[1] = 0;

				DisStruct.nodeConnectivity.insert( std::make_pair(connectivity, nodeInfo)  );
				
				isMeshed = 0; //Do not allow the remeshing to end
				break; //Break out and repeat the loop if remeshing was needed

				assert(DisStruct.nodePositions.size()==DisStruct.nodeConnectivity.size()); //Make sure we maintain the required lengths

			}

			

			//If you successfully make it out without needed to remesh, end the loop

			if(isMeshed==0)
				{
				break;	//Sanity check
				}
		}
	
	}

}

void updateNodes(DislocationStructure& DisStruct, int nodeNum1, int nodeNum2)
//Function to add or remove nodes as needed
{
	int intoMin;
	int minNode;
	//int minNode = min(nodeNum1, nodeNum2);

	int outOfMax;
	int maxNode;
	//int maxNode = max(nodeNum1, nodeNum2);

		//First, figure out orientation of the two nodes (minimum node being the one closer to node 0)
		for( auto const& [key, val] : DisStruct.nodeConnectivity)
		{
			if(key[0]==nodeNum1 && key[1]==nodeNum2) 
				{				
					minNode = nodeNum1;
					maxNode = nodeNum2;
				}	

			if(key[0]==nodeNum2 && key[1]==nodeNum1) 
				{				
					minNode = nodeNum2;
					maxNode = nodeNum1;
				}
		}

		//Find the 1st degree linkages
		for( auto const& [key, val] : DisStruct.nodeConnectivity)
		{
			if(key[0]==maxNode) //Coming out of the end
				outOfMax = key[1];

			if(key[1]==minNode) //Coming out of the end
				intoMin = key[0];	
		}

	int intoMinLink;
	int outOfMaxLink;

		//Find the 2nd degree linkages
		for( auto const& [key, val] : DisStruct.nodeConnectivity)
		{
			if(key[0]==outOfMax) //Coming out of the end
				outOfMaxLink = key[1];

			if(key[1]==intoMin) //Coming out of the end
				intoMinLink = key[0];	
		}
	/////////////////////////////////////////////////
	//// All calculations for 1st node (minNode) ////
	/////////////////////////////////////////////////

	double x_dif = abs(DisStruct.nodePositions[minNode][0] - DisStruct.nodePositions[intoMin][0]);  
	double z_dif = abs(DisStruct.nodePositions[minNode][2] - DisStruct.nodePositions[intoMin][2]);

	//Add a node if there is curvature
	if (x_dif>0 && z_dif>0)
	{
		int newNodeNum = DisStruct.nodePositions.size(); //ID number of new node to add

		//Create and insert the new node
		std::array<double,6> position;
		position[0]= DisStruct.nodePositions[intoMin][0];
		position[1]= DisStruct.nodePositions[intoMin][1];
		position[2]= DisStruct.nodePositions[minNode][2];
		position[3] = 0;
		position[4] = 0;
		position[5] = 0;

		DisStruct.nodePositions.insert( std::make_pair(newNodeNum, position)  );

		//Update the node connectivity 

		//Delete old connectivity
		std::array<int, 2> nodes;

		nodes[0] = intoMin;
		nodes[1] = minNode;

		DisStruct.nodeConnectivity.erase(nodes); //erase the old connectivity value

		//Insert the new connectivity values
		std::array<int,2> connectivity;
		connectivity[0] = intoMin;
		connectivity[1] = newNodeNum;

		std::array<int,2> nodeInfo;
		nodeInfo[0] = 0;
		nodeInfo[1] = 0;

		DisStruct.nodeConnectivity.insert( std::make_pair(connectivity, nodeInfo)  );

		//Second connectivity value
		connectivity[0] = newNodeNum;
		connectivity[1] = minNode;

		DisStruct.nodeConnectivity.insert( std::make_pair(connectivity, nodeInfo)  );

	}

	//Remove a node if two nodes now occupy the same position,
	if(x_dif==0 && z_dif==0)
	{

		//Remove the old node
		DisStruct.nodePositions.erase(intoMin);

		//Update the node connectivity//

		//Delete old connectivity
		std::array<int, 2> nodes;

		nodes[0] = intoMin;
		nodes[1] = minNode;

		DisStruct.nodeConnectivity.erase(nodes); //erase the old connectivity value

		nodes[0] = intoMinLink;
		nodes[1] = intoMin;

		DisStruct.nodeConnectivity.erase(nodes); //erase the old connectivity value

		//Insert the new connectivity values
		std::array<int,2> connectivity;
		connectivity[0] = intoMinLink;
		connectivity[1] = minNode;

		std::array<int,2> nodeInfo;
		nodeInfo[0] = 0;
		nodeInfo[1] = 0;

		DisStruct.nodeConnectivity.insert( std::make_pair(connectivity, nodeInfo)  );

	}

	/////////////////////////////////////////////////
	//// All calculations for 2nd node (maxNode) ////
	/////////////////////////////////////////////////

	x_dif = abs(DisStruct.nodePositions[maxNode][0] - DisStruct.nodePositions[outOfMax][0]);  
	z_dif = abs(DisStruct.nodePositions[maxNode][2] - DisStruct.nodePositions[outOfMax][2]);

	//Add a node if there is curvature
	if (x_dif>0 && z_dif>0)
	{
		int newNodeNum = DisStruct.nodePositions.size(); //ID number of new node to add

		//Create and insert the new node
		std::array<double,6> position;
		position[0]= DisStruct.nodePositions[outOfMax][0];
		position[1]= DisStruct.nodePositions[outOfMax][1];
		position[2]= DisStruct.nodePositions[maxNode][2];
		position[3] = 0;
		position[4] = 0;
		position[5] = 0;

		DisStruct.nodePositions.insert( std::make_pair(newNodeNum, position)  );

		//Update the node connectivity 

		//Delete old connectivity
		std::array<int, 2> nodes;

		nodes[0] = maxNode;
		nodes[1] = outOfMax;

		DisStruct.nodeConnectivity.erase(nodes); //erase the old connectivity value

		//Insert the new connectivity values
		std::array<int,2> connectivity;
		connectivity[0] = maxNode;
		connectivity[1] = newNodeNum;

		std::array<int,2> nodeInfo;
		nodeInfo[0] = 0;
		nodeInfo[1] = 0; 

		DisStruct.nodeConnectivity.insert( std::make_pair(connectivity, nodeInfo)  );


		//Second connectivity value
		connectivity[0] = newNodeNum;
		connectivity[1] = outOfMax;

		DisStruct.nodeConnectivity.insert( std::make_pair(connectivity, nodeInfo)  );

	}

	//Remove a node if two nodes now occupy the same position,
	if(x_dif==0 && z_dif==0)
	{

		//Remove the old node
		DisStruct.nodePositions.erase(outOfMax);

		//Update the node connectivity//

		//Delete old connectivity
		std::array<int, 2> nodes;

		nodes[0] = outOfMax;
		nodes[1] = outOfMaxLink;

		DisStruct.nodeConnectivity.erase(nodes); //erase the old connectivity value

		nodes[0] = maxNode;
		nodes[1] = outOfMax;

		DisStruct.nodeConnectivity.erase(nodes); //erase the old connectivity value

		//Insert the new connectivity values
		std::array<int,2> connectivity;
		connectivity[0] = maxNode;
		connectivity[1] = outOfMaxLink;

		std::array<int,2> nodeInfo;
		nodeInfo[0] = 0;
		nodeInfo[1] = 0;

		DisStruct.nodeConnectivity.insert( std::make_pair(connectivity, nodeInfo)  );

	}


}


void findVacancyIntersections(vector<Vacancy> &vacancyArray, DislocationStructure& DisStruct)
//Function to find all the absorptions of the vacancies to segments, then update the nodal positons
{
	int vacsAbsorbed=0;

	double distanceInfo[3];

	for(int i = 0; i<vacancyArray.size(); i++)
			{
				findClosestSegment(vacancyArray[i].position, DisStruct,  distanceInfo);

				//std::cout << "Vac " << vacancyArray[i].vacIDnum << " closest segment is " << distanceInfo[0] << " b" << std::endl; 

				//Only allow edge segments to absorb 
				if (distanceInfo[0]<distToAbsorbption && abs(DisStruct.nodePositions[distanceInfo[1]][2]-DisStruct.nodePositions[distanceInfo[2]][2])>0  ) 
					{
					/////////////////////////////////////////////////////////////////////////////////////////
					//Execute vacancy absorption and place the vacancy in a new position	
					/////////////////////////////////////////////////////////////////////////////////////////				
					
					vacsAbsorbed++;
					RunningVacAbsorbed++;
					
					if(useEmission==0)
					{
						//Absorb the vacancy and randomly place it somewhere else
						double pos[3];

						std::uniform_real_distribution<double> L1Dist(-L1/2,L1/2);
						std::uniform_real_distribution<double> L2Dist(-L3/2,L2/2);
						std::uniform_real_distribution<double> L3Dist(-L3/2,L3/2);

						pos[0] = L1Dist(generator);	
						pos[1] = L2Dist(generator);	
						pos[2] = L3Dist(generator);	

						for(int j = 0; j<3;j++)
							vacancyArray[i].position[j] = pos[j]*b;

						vacancyArray[i].initializeVacancy(DisStruct);

						RunningVacIDnum++;
						vacancyArray[i].vacIDnum=RunningVacIDnum;
					}
					else
					{
						vacancyArray.erase(vacancyArray.begin()+i); //Absorb and remove the vacancy if vacancy emission is turned on
					}

					/////////////////////////////////////////////////////////////////////////////////////////				
					//Execute CLIMB of the segment that the vacancy was absorbed into - via movement of the two segment
					/////////////////////////////////////////////////////////////////////////////////////////

					double segmentLength = L1/(numNodes-1); //Individual dislocation segment length in b units

					double atomicVolumeinb = atomicVolume*pow(b,-3); //Atomic volume in b units

					double h = ( atomicVolumeinb*(1+volumetricStrain)/ segmentLength); //Climb height in b units according to volume swept out by absorption of one vacancy

					//DisStruct.nodePositions[distanceInfo[1]][0]+=3; //Incorrect climb amount - used for visalization purposes
					DisStruct.nodePositions[distanceInfo[1]][0]+=h; //Correct Climb amount


					//DisStruct.nodePositions[distanceInfo[2]][0]+=3; //Incorrect climb amount - used for visalization purposes
					DisStruct.nodePositions[distanceInfo[2]][0]+=h; //Correct Climb amount
				
					
					remesh(DisStruct); //Function to remesh the dislocation structure

					}
			}

		std::cout << vacsAbsorbed << " vacancies absorbed" << std::endl;
		std::cout << RunningVacAbsorbed << " total vacancies absorbed" << std::endl;
}






#endif


