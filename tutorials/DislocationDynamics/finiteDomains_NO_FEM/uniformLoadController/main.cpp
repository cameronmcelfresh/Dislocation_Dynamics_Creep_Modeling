/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

// Define the non-singluar method used for calculations
#define _MODEL_NON_SINGULAR_DD_ 1 // 0 classical theory, 1 Cai's regularization method, 2 Lazar's regularization method

// Select the external load controller (if nothing is defined DummyExternalLoadController.h is used)
//#define ExternalLoadControllerFile <model/DislocationDynamics/ExternalLoadControllers/UniformExternalLoadController.h>
//#define ExternalLoadControllerFile <model/DislocationDynamics/ExternalLoadControllers/SampleUserStressController.h>
//#define ExternalLoadControllerFile <model/DislocationDynamics/ExternalLoadControllers/SequentialTorsionTensionController.h>
//#define ExternalLoadControllerFile <model/DislocationDynamics/ExternalLoadControllers/ClockIndentationController.h>

//#include <model/DislocationDynamics/DislocationNetwork.h>
#include <DefectiveCrystal.h>

using namespace model;

int main (int argc, char* argv[])
{
    
	//Noraml DD run
	if(useParametricStudy ==0)
	{
	
		//Create the DislocationNetwork object
		//DislocationNetwork<3,0,Hermite> DN(argc,argv);
		DefectiveCrystal<3,0,Hermite> DC(argc,argv);
		    
		// Run time steps
		DC.runGlideSteps();

	}

	//Loop through the parameters if a parametric study is chosen
	if(useParametricStudy==1)
	{

		for(int iTemp =0; iTemp<trials+1; iTemp++) //Temperature variation loop
		{

			Temp = minTemp + ((maxTemp-minTemp)/trials)*iTemp; //Temperature to use during this trial
			rewriteTemperature(Temp); //Update the W material file with the new temperature

			for(int iVacConc= 0; iVacConc<trials+1; iVacConc++) //Vacancy concentration variation loop
			{


				double thermalVacancyConcentration = (density*Na*1*exp(-1*vacFormationEnergy/(Temp*Kb))/molarMass )  * pow(100,3); //Stress-independent thermal vacancy concentration [vac/m^3]

				vacancyConcentration = thermalVacancyConcentration + ((maxVacConcentration*thermalVacancyConcentration - thermalVacancyConcentration)/trials)*iVacConc; //Vacancy concentration to use during this trial
				

				if(useStress==1) //Only enter the stress variation loops if there is applied stress
				{
					for(int iStress = 0; iStress<trials+1; iStress++) //sigma_xx stress variation loop
					{

						appliedPressure = minPressure + ((maxPressure-minPressure)/trials)*iStress; //sigma_xx to use during this trial

						for(int repeatCycles=0; repeatCycles<repeatTrials; repeatCycles++) //Repeat the simulations however many times desired
						{	


							resetRunningVariables(); //Reset all of the counters prior to running each new trial

							//////Reset the counters each trial!!//////
							RunningVacAbsorbed=0; //Running counter of how many vacancies have been absorbed
							RunningLastVacNumber = 0; //Running counter of how many vacancies were previously absorbed before last calculation of dislocation velocity
							RunningVacEmitted=0; //Running counter of how many vacancies have been emitted
							RunningLastVacEmitted=0; //Running counter of how many vacancies were previously emitted before last calculation of dislocation velocity
							RunningVacIDnum = 0; //Running counter for the ID of vacancies
							totalGlobalTime = 0; //Running counter for the total gloabl time
							lastTotalGlobalTime = 0; //Running counter for the previous gloabl time to be used in velocity calculation
							lastDistanceMoved = 0; //Placeholder for the last distance moved
							//////////////////////////////////////////
								
							//Create the DislocationNetwork object
							DefectiveCrystal<3,0,Hermite> DC(argc,argv);
										    
							// Run time steps
							DC.runGlideSteps();

						}
		
					}
				}
				else //Skip the stress variations if there is no applied stress
				{
					for(int repeatCycles=0; repeatCycles<repeatTrials; repeatCycles++) //Repeat the simulations however many times desired
						{	

						resetRunningVariables(); //Reset all of the counters prior to running each new trial

						//////Reset the counters each trial!!//////
						RunningVacAbsorbed=0; //Running counter of how many vacancies have been absorbed
						RunningLastVacNumber = 0; //Running counter of how many vacancies were previously absorbed before last calculation of dislocation velocity
						RunningVacEmitted=0; //Running counter of how many vacancies have been emitted
						RunningLastVacEmitted=0; //Running counter of how many vacancies were previously emitted before last calculation of dislocation velocity
						RunningVacIDnum = 0; //Running counter for the ID of vacancies
						totalGlobalTime = 0; //Running counter for the total gloabl time
						lastTotalGlobalTime = 0; //Running counter for the previous gloabl time to be used in velocity calculation
						lastDistanceMoved = 0; //Placeholder for the last distance moved
						//////////////////////////////////////////
								
						//Create the DislocationNetwork object
						DefectiveCrystal<3,0,Hermite> DC(argc,argv);
										    
						// Run time steps
						DC.runGlideSteps();

						}

				}
				


			}//end vacancy concetration variation
		}//end temperature variation
	}//end parametric study
    
    return 0;
}
