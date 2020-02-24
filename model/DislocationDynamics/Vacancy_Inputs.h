#include <DefectiveCrystal.h>
#include <DislocationNetwork.h>
#include <cmath>
#include <random>
#include <string>


#ifndef Vacancy_Inputs_H_
#define Vacancy_Inputs_H_

using namespace std;

////////Physical Constants////////
double Kb = 1.3806*pow(10,-23); //Boltzmann Constant
double Na = 6.022*pow(10,23); //Avogadro's number in atoms/mole

////////Material Constants////////
string materialName("Fe"); //Material name -- must match with a material file listed in ../../MaterialsLibrary/
double density = 7.86;  //Bulk density of material [g/cm^3]
double molarMass = 55.845; //Molar mass [g/mol]
double vacFormationEnergy = 2.7237*pow(10,-19); //Energy per vacancy formation [Joules]

double b = 0.25*pow(10,-9); //Burger's vector [m]
double G = 80*pow(10,9); //Shear modulus in [Pa]
double poissonRatio = 0.29; //Poisson's Ratio
double Eo = 9.6*pow(10,-20); //Vacancy Migration Energy [Joules]
double v = 1*pow(10,12); //Vibrational frequency [Hz]

double atomicVolume = 0.77*pow(b,3); //Atomic volume [m^3]
double volumetricStrain = 0.1; //Volumetric strain
double nearestNeighbors = 8; //8 Nearest neighbors for bcc
double correlationFactor = 0.78; //Correlation factor 

////////Simulation Constants////////
int numNodes = 14; //Number of nodes to discretize across the straight dislocation that runs across the simulation box
double fractZOverhang = .1; //Fraction of the total z-length that the z-length straight dislocaiton extends above and below the simulation box
double fractXOverhang = .2; //Fraction of the total x-length that the dislocation loop extends beyond the simulation box
int DisStructUpdateFreq = 5; //Number of runID trials between updating the DislocationStructure and absorbing vacancies
int PrintStatsFreq = 1000; //Frequency of outputting the test statistics
int PrintEVLFreq = 1000; //Frequency of outputting the evl file
int PrintCPUTime = 1; //Whether or not to print the running CPU time for the simulation in results/stats.txt

int kmcAcceptableError = 30; //Acceptable percentage error for self-guess of kMC timestep
int outputGlobalTimeStep = 0; //1=true or 0=false for printing globaltimestep file -- printed in "test" folder
int outputV = 0; //1 or 0 for printing V_#.txt files -- printed in "evl1" folder
int outputEVL=0; //1 or 0 for printing evl_#.txt files -- printed in "evl1" folder

////Do not change//////
double totalGlobalTime = 0; //Running counter for the total gloabl time
double lastTotalGlobalTime = 0; //Running counter for the previous gloabl time to be used in velocity calculation
double lastDistanceMoved = 0;

////////Mesh-Related////////
int meshType = 0; // 0 = rectangular prism mesh, 1 = cylindrical mesh

double L1=1632; //the side length of the cube, in units of Burgers vector
double L2=1510; //the side length of the cube, in units of Burgers vector
double L3=1520; //the side length of the cube, in units of Burgers vector

//If meshType==0
//L1 = x-length
//L2 = y-length
//L3 = z-length

//If meshType==0
//L1 = radius
//L2 = z-length

////////Physically Relevant Simulation Constants////////
double Temp = 800; //Temperature (Kelvin)
double minVacJump = b/100; //Minimum jump distance of vacancy for MC purposes
double maxVacJump = b; //Maximum jump distance of vacancy for MC purposes

////////Vacancy Related////////
double distToAbsorbption = 2; //Distance in Burger's units for a vacancy to be absorbed by a dislocation
double vacancyConcentration =0; //2e21; //Volumetric Concentration of vacancies in vac/m^3 - overrides the calculation of the number of vacancies if not vacancyConcentration=0
int useEmission = 1; //0=do not use emission of vacancaies -- randomly replace them each time they are absorbed, 1=do not automatically replace vacancies -- have segments emit vacancies and negatively climb
int useDiscreteEmission = 1; //1 = use discrete emission from local concentration of vacancies around a dislocation segment, 0=use a global average of vacancy concentration for vacancy emission

////DO NOT CHANGE////
int RunningVacAbsorbed=0; //Running counter of how many vacancies have been absorbed
int RunningLastVacNumber = 0; //Running counter of how many vacancies were previously absorbed before last calculation of dislocation velocity
int RunningVacIDnum = 0; //Running counter for the ID of vacancies
int RunningVacEmitted=0; //Running counter of how many vacancies have been emitted
int RunningLastVacEmitted=0; //Running counter of how many vacancies were previously emitted before last calculation of dislocation velocity
int RunningBalancedVacs = 0; //Running counter of how many vancacies were added (+) or subtracted (-) from the system to keep the concentration constant
double RunningCPUTime = 0; //Running counter of the total CPU time elapsed
int vacNum; //Number of vacancies - to be calculated later

////////Random Generators////////
random_device generator; //Generator to select random numbers throughout the simulation

////////For Vacancy Movement////////
uniform_real_distribution<double> ZeroOnedistribution(0, 1); //Uniform distribution to be used throughout program
uniform_real_distribution<double> Stepdistribution(minVacJump, maxVacJump); //Uniform distribution for step distances to be used throughout the program

////////For Random Vacancy Placement In Mesh//////// -- had to explicitly redefine each generator within the vacancies/defectiveCrystal file due to dynamicBoxResizing 
//std::uniform_real_distribution<double> L1Dist(-L1/2,L1/2);
//std::uniform_real_distribution<double> L2Dist(-L3/2,L2/2);
//std::uniform_real_distribution<double> L3Dist(-L3/2,L3/2);

////////Parametric Testing of Climb Parameters////////
int useParametricStudy = 1; //0=run non-parametric study, 1 = run parametric study by varying the variables listed below
int DDSimulationText=0; //0=minimal output from the DD module, 1=full DD and vacancy simulation text output

int dynamicBoxResizing = 1; //0=do not dynamically resize the box to maintain a constant vacancy number, 1 = resize box as to keep a constant vacancy number throughout the parametric study
int constVacNumber = 100; //Number of vacancies to consider for each trial of the parametric stdy IF dynamicBoxResizing==1

int requireConstantVacs = 1; //0=do not add/remove vacancies to keep the pre-set vacancy concentration, 1=ensure that the pre-set vacancy concentration is kept by adding/removing vacancies

int tempTrials = 4; //Number of intermediary temperature trials to do within min and max bounds
int concTrials = 4; //Number of intermediary concentration multiplier trials to do within min and max bounds
int stressTrials = 4; //Number of intermediary stress trials to do within min and max bounds
int trials = 4; //Number of trials to do within the min and max bounds. For example, if minTemp=100, maxTemp=200, and trials = 3, simulations would be run at 100, 150, and 300.
int repeatTrials = 3; //How many times to repeat each simulation

double minTemp = 800; //Minimum temperature to use in parametric study [Kelvins]
double maxTemp = 1500; //Maximum temperature to use in parametric study [Kelvins]

double minVacConcentration = 1; //Minimum vacancy concentration to use in parametric study  - must be a multiplier of the thermal vacancy concentration[vacs/m^3] 
double maxVacConcentration = 100; //Maximum vacancy concentration to use in parametric study - must be a multiplier of the thermal vacancy concentration [vacs/m^3]
//If the thermal vacancy concentration is 1E10 and maxVacConcentration=5, the maximum vacancy concentration will be 5E10.

int useStress = 1;  //0 = do not use applied stress in the parametric study, 1 = use applied stress in the parametric study
double minPressure = -100*pow(10,6); //Minimum stress [pressure] to use in parametric study [Pa]
double maxPressure = 100*pow(10,6);  //Maximum stress [pressure] to use in parametric study [Pa]
double appliedPressure = 0; //DO NOT CHANGE -- this is the default (thermal) pressure used to calculate the enthalpy of vacancy formation




#endif

