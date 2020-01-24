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
    
	//std::cout.rdbuf(nullptr); //Supress output
	std::streambuf* cout_sbuf = std::cout.rdbuf(); // save original sbuf
	std::ofstream   fout("/dev/null");

	if(DDSimulationText==0)
	{
	std::cout.rdbuf(fout.rdbuf()); // redirect 'cout' to a 'fout' -- mute all DD console output text
	}

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

		double trialCounter =0; //Counter to estimate time remaining
		auto runningSeconds= std::chrono::system_clock::now(); //Clock to count length of trials

		for(int iTemp =0; iTemp<trials+1; iTemp++) //Temperature variation loop
		{

			Temp = minTemp + ((maxTemp-minTemp)/trials)*iTemp; //Temperature to use during this trial
			rewriteTemperature(Temp); //Update the material file with the new temperature

			//Set the number of DD steps to run (Nsteps variable) based on the simulation temperature
			if(Temp==800)
				{
				double newNsteps = 10000;
				rewriteDD(newNsteps);
				}

			if(Temp>1200)
				{
				double newNsteps = 1500100;
				rewriteDD(newNsteps);
				}

				
			for(int iVacConc= 0; iVacConc<trials+1; iVacConc++) //Vacancy concentration variation loop
			{
				

				if(useStress==1) //Only enter the stress variation loops if there is applied stress
				{
					for(int iStress = 0; iStress<trials+1; iStress++) //sigma_xx stress variation loop
					{

						int TotalTrials = pow(trials,3)*repeatTrials; //Total number of trials that will be run

						appliedPressure = minPressure + ((maxPressure-minPressure)/trials)*iStress; //sigma_xx to use during this trial


						//Once the pressure has been defined, calculate the corresponding temperature-pressure dependent vacancy concentration
						double formationEnthalpy = vacFormationEnergy - appliedPressure*(1+volumetricStrain)*atomicVolume; //Stress-dependent vacancy formation energy
						double thermalVacancyConcentration = (density*Na*1*exp(-1*formationEnthalpy/(Temp*Kb))/molarMass )  * pow(100,3); //Stress-independent thermal vacancy concentration [vac/m^3]
						vacancyConcentration = thermalVacancyConcentration + ((maxVacConcentration*thermalVacancyConcentration - thermalVacancyConcentration)/trials)*iVacConc; //Vacancy concentration to use during this trial

						for(int repeatCycles=0; repeatCycles<repeatTrials; repeatCycles++) //Repeat the simulations however many times desired
						{	

							resetRunningVariables(); //Reset all of the counters prior to running each new trial
								
							//Create the DislocationNetwork object
							DefectiveCrystal<3,0,Hermite> DC(argc,argv);
										    
							// Run time steps
							DC.runGlideSteps();

							double secs = std::chrono::duration<double>(std::chrono::system_clock::now()-runningSeconds).count(); //Total seconds elapsed
							double hours = secs/(60*60); //Hours Elapsed

							trialCounter++;
							double percentage_elapsed = trialCounter/(pow(trials,3)*repeatTrials); //Percentage of the trial that has elapsed

							double hours_left = hours/percentage_elapsed-hours; //Total hours left in the simulation

							if(DDSimulationText==0)
								{
   								std::cout.rdbuf(cout_sbuf); // restore the original stream buffer
								std::cout<<"Completed T = " << Temp << " [K] , Vac Concentration = " << vacancyConcentration << " [vacs/m^3] , " << "Pressure = " << appliedPressure << " [Pa] || "; 

								if(hours_left<1)
									std::cout<< percentage_elapsed*100 <<"%" << " done ---> ~" << hours_left*60 << " minutes remain" << std::endl;
								else
									std::cout<< percentage_elapsed*100 <<"%" << " done ---> ~" << hours_left << " hours remain" << std::endl;

								std::cout.rdbuf(fout.rdbuf()); // redirect 'cout' to a 'fout'
								}

						}
		
					}
				}
				else //Skip the stress variations if there is no applied stress
				{
					for(int repeatCycles=0; repeatCycles<repeatTrials; repeatCycles++) //Repeat the simulations however many times desired
						{	


						resetRunningVariables(); //Reset all of the counters prior to running each new trial
								
						//Create the DislocationNetwork object
						DefectiveCrystal<3,0,Hermite> DC(argc,argv);
										    
						// Run time steps
						DC.runGlideSteps();

						double secs = std::chrono::duration<double>(std::chrono::system_clock::now()-runningSeconds).count(); //Total seconds elapsed
						double hours = secs/(60*60); //Hours Elapsed

						trialCounter++;
						double percentage_elapsed = trialCounter/(pow(trials,2)*repeatTrials);//Percentage of the trial that has elapsed

						double hours_left = hours/percentage_elapsed-hours; //Total hours left in the simulation

						if(DDSimulationText==0)
							{
							std::cout.rdbuf(cout_sbuf); // restore the original stream buffer
							std::cout<<"Completed T = " << Temp << " [K] , Vac Concentration = " << vacancyConcentration << " [vacs/m^3] || ";

							if(hours_left<1)
								std::cout<< percentage_elapsed*100 <<"%" << " done ---> ~" << hours_left*60 << " minutes remain" << std::endl;
							else
								std::cout<< percentage_elapsed*100 <<"%" << " done ---> ~" << hours_left << " hours remain" << std::endl;

							std::cout.rdbuf(fout.rdbuf()); // redirect 'cout' to a 'fout'
							}

						}

				}
				


			}//end vacancy concetration variation
		}//end temperature variation
	}//end parametric study
    
    return 0;
}
