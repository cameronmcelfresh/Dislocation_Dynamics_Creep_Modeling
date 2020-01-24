/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2017 by Yinan Cui <cuiyinan@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DefectiveCrystal_H_
#define model_DefectiveCrystal_H_

#include <iostream>
#include <vector>
#include <memory>

//#ifndef ExternalLoadControllerFile
//#define ExternalLoadControllerFile <model/DislocationDynamics/ExternalLoadControllers/DummyExternalLoadController.h>
//#endif
//#include ExternalLoadControllerFile

#include <DefectiveCrystalParameters.h>
#include <DislocationNetwork.h>
#include <CrackSystem.h>
#include <UniformExternalLoadController.h>

//#include <Test.h>
#include <Vacancies.h>


namespace model
{
    
    template <int _dim, short unsigned int corder, typename InterpolationType>
    class DefectiveCrystal
    {
        
    public:
        static constexpr int dim=_dim; // make dim available outside class
        typedef DefectiveCrystal<dim,corder,InterpolationType> DefectiveCrystalType;
        typedef DislocationNetwork<dim,corder,InterpolationType> DislocationNetworkType;
        typedef CrackSystem<dim> CrackSystemType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef BVPsolver<dim,2> BVPsolverType;
        typedef typename BVPsolverType::ElementType ElementType;
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Cameron's Variables
	//static constexpr int vacNum=100; //////////////////////////////////////////////// CAMERON'S CODE
	vector<Vacancy> vacancies; //all vacancies
	DislocationStructure DisStruct; //dislocation structure object for coupling dislocations and mesh with vacancies

	
	//static constexpr int DisStructUpdateFreq = 10; //Number of runID trials between updating the DislocationStructure 
	//static constexpr int vacNumToClimb = 2; //Number of vacancies that a dislocation segment needs to absorb in order to climb

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        DefectiveCrystalParameters simulationParameters;
        
        //        long int runID;
        
        
        
        const SimplicialMesh<dim> mesh;
        const std::vector<VectorDim> periodicShifts;
        const Polycrystal<dim> poly;
        const std::unique_ptr<DislocationNetworkType> DN;
        const std::unique_ptr<CrackSystemType> CS;
        const std::unique_ptr<BVPsolverType> bvpSolver;
        const std::unique_ptr<ExternalLoadControllerBase<dim>> externalLoadController;
        
        
        /**********************************************************************/
        static std::unique_ptr<ExternalLoadControllerBase<dim>> getExternalLoadController(const DefectiveCrystalParameters& params,
                                                                                          const DefectiveCrystalType& dc,
                                                                                          const long int& rID)
        {
            
            if(params.simulationType==DefectiveCrystalParameters::FINITE_FEM)
            {
                return std::unique_ptr<ExternalLoadControllerBase<dim>>(nullptr);
            }
            else
            {
                if(params.externalLoadControllerName=="UniformExternalLoadController")
                {
                    return std::unique_ptr<ExternalLoadControllerBase<dim>>(new UniformExternalLoadController<DefectiveCrystalType>(dc,rID));
                }
                else if(params.externalLoadControllerName=="None")
                {
                    return std::unique_ptr<ExternalLoadControllerBase<dim>>(nullptr);
                }
                else
                {
                    model::cout<<"Unknown externalLoadController name "<<params.externalLoadControllerName<<". Use 'None' for no-controller. EXITING."<<std::endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
        
        /**********************************************************************/
        static std::vector<VectorDim> getPeriodicShifts(const SimplicialMesh<dim>& m,
                                                        const DefectiveCrystalParameters& params)
        {
            // Set up periodic shifts
            std::vector<VectorDim> temp;
            if(params.simulationType==DefectiveCrystalParameters::PERIODIC_IMAGES)
            {
                const VectorDim meshDimensions(m.xMax()-m.xMin());
                model::cout<<"meshDimensions="<<meshDimensions.transpose()<<std::endl;
                for(int i=-params.periodicImages_x;i<=params.periodicImages_x;++i)
                {
                    for(int j=-params.periodicImages_y;j<=params.periodicImages_y;++j)
                    {
                        for(int k=-params.periodicImages_z;k<=params.periodicImages_z;++k)
                        {
                            const Eigen::Array<int,dim,1> cellID((Eigen::Array<int,dim,1>()<<i,j,k).finished());
                            temp.push_back((meshDimensions.array()*cellID.template cast<double>()).matrix());
                        }
                    }
                }
            }
            else
            {
                temp.push_back(VectorDim::Zero());
            }
            
            model::cout<<"periodic shift vectors:"<<std::endl;
            for(const auto& shift : temp)
            {
                model::cout<<shift.transpose()<<std::endl;
                
            }
            
            return temp;
            
        }
        
        
        
        /**********************************************************************/
        void updateLoadControllers(const long int& runID, const bool& isClimbStep)
        {/*! Updates bvpSolver using the stress and displacement fields of the
          *  current DD configuration.
          */
            const int quadraturePerTriangle=37;
            if(bvpSolver)
            {
                if (!(runID%bvpSolver->stepsBetweenBVPupdates))
                {// enter the if statement if use_bvp!=0 and runID is a multiple of use_bvp
                    model::cout<<"		Updating elastic bvp... "<<std::endl;
                    bvpSolver->template assembleAndSolve<DislocationNetworkType,quadraturePerTriangle>(*DN, isClimbStep);
                }
            }
            if (externalLoadController)
            {
                externalLoadController->update(runID);
            }
        }
        
    public:
        
        
        /**********************************************************************/
        DefectiveCrystal(int& argc, char* argv[]) :
        /* init */ simulationParameters(argc,argv)
        /* init */,mesh(TextFileParser("inputFiles/DD.txt").readScalar<int>("meshID",true))
        /* init */,periodicShifts(getPeriodicShifts(mesh,simulationParameters))
        /* init */,poly("./inputFiles/polycrystal.txt",mesh)
        /* init */,DN(simulationParameters.useDislocations? new DislocationNetworkType(argc,argv,simulationParameters,mesh,poly,bvpSolver,externalLoadController,periodicShifts,simulationParameters.runID) : nullptr)
        /* init */,CS(simulationParameters.useCracks? new CrackSystemType() : nullptr)
        //        /* init */,DN(argc,argv,simulationParameters,mesh,poly,bvpSolver,externalLoadController,periodicShifts,simulationParameters.runID)
        /* init */,bvpSolver(simulationParameters.simulationType==DefectiveCrystalParameters::FINITE_FEM? new BVPsolverType(mesh,*DN) : nullptr)
        /* init */,externalLoadController(getExternalLoadController(simulationParameters,*this,simulationParameters.runID))
        {
            assert(mesh.simplices().size() && "MESH IS EMPTY.");
            
            
            if(   simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_IMAGES
               || simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_FEM)
            {
                assert(poly.grains().size()==1 && "ONLY SINGLE-CRYSTAL PERIODIC SIMULATIONS SUPPORTED.");
                
                for(const auto& rIter : mesh.regions())
                {
                    for(const auto& pair : rIter.second->parallelFaces())
                    {
                        model::cout<<"Checking if parallel faces "<<pair.first<<"<->"<<pair.second<<" are commensurate"<<std::endl;
                        const PlanarMeshFace<dim>& face1(*rIter.second->faces().at(pair.first));
                        const PlanarMeshFace<dim>& face2(*rIter.second->faces().at(pair.second));
                        const VectorDim cc(face1.center()-face2.center());
                        const VectorDim ccc(cc.dot(face1.outNormal())*face1.outNormal());
                        
                        const LatticeDirection<dim> ld(poly.grains().begin()->second.latticeDirection(face1.outNormal()));
                        const double normRatio(ccc.norm()/ld.cartesian().norm());
                        if(std::fabs(std::round(normRatio)-normRatio)>FLT_EPSILON)
                        {
//                            std::cout<<"Face outNormal="<<std::setprecision(15)<<std::scientific<<face1.outNormal().transpose()<<std::endl;
                            std::cout<<"Mesh in direction "<< std::setprecision(15)<<std::scientific<<ld.cartesian().normalized().transpose()<<" is not commensurate for periodicity"<<std::endl;
                            std::cout<<"Mesh size in that direction must be a multiple of "<< std::setprecision(15)<<std::scientific<<ld.cartesian().norm()<<std::endl;
                            std::cout<<"Size detected="<< std::setprecision(15)<<std::scientific<<ccc.norm()<<std::endl;
                            std::cout<<"Closest commensurate size="<< std::setprecision(15)<<std::scientific<<std::round(normRatio)*ld.cartesian().norm()<<std::endl;
                            assert(false && "MESH NOT COMMENSURATE");
                        }
//                        LatticeVector<dim> lv(ccc,poly.grains().begin()->second);
                    }
                }
            }
            
            //            switch (simulationParameters.simulationType)
            //            {// Initilization based on type of simulation
            //
            //
            //                case DefectiveCrystalParameters::FINITE_NO_FEM:
            //                {
            //                    //                    externalLoadController->init(DN,runID);  // have to initialize it after mesh!
            //                    break;
            //                }
            //
            //                case DefectiveCrystalParameters::FINITE_FEM:
            //                {
            //                    //                    bvpSolver->init(DN);
            //                    break;
            //                }
            //
            //                case DefectiveCrystalParameters::PERIODIC_WITH_IMAGES:
            //                {
            //
            ////
            ////                    const VectorDim meshSize(mesh.xMax()-mesh.xMin());
            ////                    for(int d=0;d<dim;++d)
            ////                    {
            ////                        VectorDim v(VectorDim::Zero());
            ////                        v(d)=1*meshSize(d);
            ////                        LatticeVector<dim> lv(v,poly.grains().begin()->second);
            ////                    }
            //                    break;
            //                }
            //
            //                default:
            //                {
            //                    model::cout<<"simulationType MUST BE 0,1, or 2. EXITING."<<std::endl;
            //                    exit(EXIT_FAILURE);
            //                    break;
            //                }
            //            }
       

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//CAMERON's part of the constructor
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	

	//////Calculate the number of vacancies needed//////
	if (vacancyConcentration==0) //Calculate the vacancy concentration from temperature and activation energy if not pre-set
		{
		double formationEnthalpy = vacFormationEnergy; //Stress-independent vacancy formation energy

		
		vacancyConcentration = ( density*Na*1*exp(-1*formationEnthalpy/(Temp*Kb))/molarMass ) * pow(100,3); //Temperature-dependednt vacancy concentration [vacs/m^3]
		std::cout << "Vacancy Concentration = " << vacancyConcentration << std::endl;		
		}

	if (meshType==0)
		{
		vacNum = round(vacancyConcentration*pow(b,3)*L1*L2*L3); //Rectangular simulation mesh
		std::cout << "vacNum = " << vacNum << std::endl;		
		}
	if (meshType==1)
		vacNum = round(vacancyConcentration*3.1415*pow(L1*b,2)*L2*b); //Cylindrical simulation mesh


	if( (dynamicBoxResizing==1) && (useParametricStudy==1) ) //If doing a parametric study with dynamic resizing, scale the side lengths to achieve a vac concentration without changing the number of vacancies
		{
		vacNum = constVacNumber; //Constant number of vacancies for each trial
		
		//Resize L1, L2, and L3 but keep the original ratio constant
		double newVolume = vacNum * pow(vacancyConcentration, -1) *  pow(b,-3) ; //Volume required of the new simulation mesh, in b units

		double sideScalingFactor = pow(newVolume/(L1*L2*L3), 0.3333); //Factor to multiply each side length by

		//Resize L1, L2, L3
		L1 = L1*sideScalingFactor;
		L2 = L2*sideScalingFactor;
		L3 = L3*sideScalingFactor;

		}

	std::cout << "Test 1 - resizing vacancies " << std::endl;
	std::cout << vacNum << " vacancies" << std::endl;
	vacancies.resize(vacNum); //Resize the vacancies vector to the correct size
	std::cout << "Test 2 - resizing vacancies" << std::endl;
	////////////////////////////////////////////////

	////Updating the physical nodes of the dislocation network after absorbing vacancies - parametric study////
	if (useParametricStudy==1) 
		{
		DisStruct.createNodesAndSegments(numNodes,fractZOverhang, fractXOverhang); //Create the dislocation network AFTER L1, L2 and L3 have been resized (if applicable)
		}


	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


 }
        
        //        /**********************************************************************/
        //        double compute_dt() const
        //        {
        //            if(DN)
        //            {
        //            switch (simulationParameters.timeIntegrationMethod)
        //            {
        //                case 0:
        //                    return DDtimeIntegrator<0>::get_dt(*DN);
        //                    break;
        //
        //                default:
        //                    assert(0 && "time integration method not implemented");
        //                    return 0;
        //                    break;
        //            }
        //            }
        //            else
        //            {
        //                assert(0 && "dt calculation not implemented");
        //                return 0;
        //            }
        //        }
        
        /**********************************************************************/
        void singleGlideStep()
        {

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////// Cameron test space


	model::cout<< "Cameron's vacancy code made it!" << std::endl;

	double pOne[3] = {0, 0, -10};
	double pTwo[3] = {0, 0, 10};
	double vacPoint[3] = {0, 4*b,0};
	double vacPoint1[3] = {0, -4,0};
	double vacPoint2[3] = {4, 0,0};
	double vacPoint3[3] = {3, 4,0};

	//std::cout << "Test 1,2,3,4 = " << returnDistanceToSegment(pOne, pTwo, vacPoint) << " , " << returnDistanceToSegment(pOne, pTwo, vacPoint1) << " , " <<returnDistanceToSegment(pOne, pTwo, vacPoint2) << " , " <<returnDistanceToSegment(pOne, pTwo, vacPoint3) << std::endl;

	double dInfo[3];

	findClosestSegment(vacPoint, DisStruct, dInfo);

	//std::cout << "The closest segment is " << dInfo[0] << " b, between nodes " << dInfo[1] << " and " << dInfo[2] << std::endl;



	//Dislocation arm test space
	//Eigen::Matrix<double,dim,1> p0(0,0,-500);
	//Eigen::Matrix<double,dim,1> p1(0,0,-250);
	//Eigen::Matrix<double,dim,1> t(0,0,1);
	//double length = 200;
	//Eigen::Matrix<double,dim,1> burgers(1,0,0);
	//Eigen::Matrix<double,dim,1> position(3,3,0);

	//StressStraight s0 = StressStraight(p0, p1, burgers);

	//model::cout << "Test 0" << std::endl;

	//MatrixDim M0 = s0.stress(position);

	//std::cout << M0 << std::endl;
	//
	
	//VectorDim positionVec(50,50,50);
	//MatrixDim stressMat(MatrixDim::Zero());
	//stress(positionVec) = stressMat; 
	//stressMat+=DN->stress(positionVec);
	//std::cout << stressMat << std::endl;
	//model::cout << DN->nodes().size() << std::endl;
	
	////////////////////////////////////////////////////////
	//////////Gradient printing for simulation box//////////
	////////////////////////////////////////////////////////

	double posGrad[3];

	//double disPosY = 1.999697977195557e+03, disPosZ= 1.192085008154480e+03, disPosX = 2.0000000000+03;
	double tempPosX, tempPosY, tempPosZ;
	double tempGradientMatrix[3][3][3];
	double xGrad;
	double yGrad;
	double zGrad;

	int numGridPoints = 12;

	double Xlength = L1*2;
	double Ylength = L2*2;
	double Zlength = L3*2;

	double gradientStep = pow(10,-9);

	double Xincrement = Xlength/numGridPoints;
	double Yincrement = Ylength/numGridPoints;
	double Zincrement = Zlength/numGridPoints;

	//double disPosX = Xlength/2, disPosY = Ylength/2, disPosZ = Zlength/2;

	if (simulationParameters.runID==10)
	{
	
	ofstream OutGradPosStream;
	OutGradPosStream.open("./test/xyz_data.txt");
	
	ofstream OutGradStream;
	OutGradStream.open("./test/gradient_data.txt");

	for(int i=0;i<numGridPoints; i++)
		{
		for(int j = 0;j <numGridPoints; j++)
			for(int k = 0; k<numGridPoints;k++)
			{
			tempPosX = - Xlength/2 + Xincrement*k;
			tempPosY = - Ylength/2 + Yincrement*i;
			tempPosZ = - Zlength/2 + Zincrement*j;

			posGrad[0] = tempPosX*b;
			posGrad[1] = tempPosY*b;
			posGrad[2] = tempPosZ*b;

			int placeHolder =1;

			calc3x3x3Gradient(posGrad,tempGradientMatrix, gradientStep, DisStruct, 1);

			xGrad = tempGradientMatrix[0][0][0] + tempGradientMatrix[0][1][1] + tempGradientMatrix[0][2][2];

			yGrad = tempGradientMatrix[1][0][0] + tempGradientMatrix[1][1][1] + tempGradientMatrix[1][2][2];

			zGrad = tempGradientMatrix[2][0][0] + tempGradientMatrix[2][1][1] + tempGradientMatrix[2][2][2];

			OutGradPosStream << tempPosX << " "  << tempPosY << " " << tempPosZ << "\n";

			OutGradStream << xGrad/pow(pow(xGrad,2)+pow(yGrad,2)+pow(zGrad,2),0.5) <<" " << yGrad/pow(pow(xGrad,2) + pow(yGrad,2)+pow(zGrad,2),0.5) << " " << zGrad/pow(pow(xGrad,2) + pow(yGrad,2)+pow(zGrad,2),0.5) << "\n";

			}
		}

/*
			for(int i = 1; i<5;i++)
			{
			posGrad[0] = 20*i;
			posGrad[1] = 20*i;
			posGrad[2] = 20*i;
			calc3x3x3Gradient(posGrad,tempGradientMatrix, gradientStep, DN);
			yGrad = tempGradientMatrix[1][0][0] + tempGradientMatrix[1][1][1] + tempGradientMatrix[1][2][2];
			zGrad = tempGradientMatrix[2][0][0] + tempGradientMatrix[2][1][1] + tempGradientMatrix[2][2][2];
			OutGradPosStream << tempPosY << " , " << tempPosZ << "\n";
			OutGradStream << yGrad/pow(pow(yGrad,2)+pow(zGrad,2),0.5) << " , " << zGrad/pow(pow(yGrad,2)+pow(zGrad,2),0.5) << "\n";
			}
*/

	OutGradPosStream.close();
	OutGradStream.close();

	}		
 
	////////////////////////////////////////////////////
	//////////End gradient printing/////////////////////
	////////////////////////////////////////////////////

	//printstress(pos, DN);

	auto tV= std::chrono::system_clock::now();

	//Initialize the vacancies to the dislocation network on runID = 1...
	if (simulationParameters.runID==1)
	{	
		std::cout << "Initializing vacancies..." << std::endl;
		std::cout << "Vacancy Concentration = " <<vacancyConcentration<< " vac/m^3" << std::endl;
		std::cout << "Number of vacancies in mesh = " <<vacNum<<  std::endl;

		//std::default_random_engine generator;
		std::uniform_real_distribution<double> L1Dist(-L1/2,L1/2);
		std::uniform_real_distribution<double> L2Dist(-L3/2,L2/2);
		std::uniform_real_distribution<double> L3Dist(-L3/2,L3/2);
	
		double pos[3];

		for(int i = 0; i<vacNum;i++)
		{

			pos[0] = L1Dist(generator);	
			pos[1] = L2Dist(generator);	
			pos[2] = L3Dist(generator);	

			for(int j = 0; j<3;j++)
				vacancies[i].position[j] = pos[j]*b;


			vacancies[i].initializeVacancy(DisStruct);
		}
		
	}

	////Populate the Dislocation Structure////
	if (simulationParameters.runID==1 && useParametricStudy==0) 
		{
		DisStruct.readNodesAndSegments(simulationParameters.runID-1);
		}


	////Updating the physical nodes of the dislocation network after absorbing vacancies - parametric study////
	if ((simulationParameters.runID-1)%DisStructUpdateFreq==0 && simulationParameters.runID>1 && useParametricStudy==1) 
		{
		findVacancyIntersections(vacancies, DisStruct);
		}

	////Updating the physical nodes of the dislocation network after absorbing vacancies - parametric study////
	if ((simulationParameters.runID-1)%DisStructUpdateFreq==0 && simulationParameters.runID>1 && useParametricStudy==0) 
		{
		findVacancyIntersections(vacancies, DisStruct, DN);
		}
	
	////Print test stats - dislocation velocity and vacancy statistics////
	if ((simulationParameters.runID-1)%PrintStatsFreq==0 && simulationParameters.runID>1) 
		{
		printStatistics(DisStruct, simulationParameters.runID);
		}


	////Print the FIRST evl file so there is a evl_0.txt////
	if ((simulationParameters.runID)==1) 
		{
		DisStruct.printevl(simulationParameters.runID-1);
		}

	////Print evl files -- dislocation nodes and segments////
	if ((simulationParameters.runID)%PrintEVLFreq==0 && simulationParameters.runID>1) 
		{
		DisStruct.printevl(simulationParameters.runID);
		}

	////For runID>1 allow the vacancies to move! - parametric study////
	if (simulationParameters.runID>1 && useParametricStudy==1)
		{
		std::cout << "Updating vacancy positions..." << std::endl;
		singleVacancyEvent(vacancies,simulationParameters.runID,DisStruct);
		}

	////For runID>1 allow the vacancies to move! - regular DD////
	if (simulationParameters.runID>1 && useParametricStudy==0)
		{
		std::cout << "Updating vacancy positions..." << std::endl;
		singleVacancyEvent(vacancies, DN, simulationParameters.runID,DisStruct);
		}

	std::cout << "Vacany updates [" << (std::chrono::duration<double>(std::chrono::system_clock::now()-tV)).count() << " sec]" << std::endl << std::endl;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////END CAMERON'S CODE///////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            model::cout<<blueBoldColor<< "runID="<<simulationParameters.runID<<" (of "<<simulationParameters.Nsteps<<")"
            /*                    */<< ", time="<<simulationParameters.totalTime;
            if(DN)
            {
                model::cout<< ": nodes="<<DN->nodes().size()
                /*                    */<< ", segments="<<DN->links().size()
                /*                    */<< ", loopSegments="<<DN->loopLinks().size()
                /*                    */<< ", loops="<<DN->loops().size()
                /*                    */<< ", components="<<DN->components().size();
            }
            model::cout<< defaultColor<<std::endl;
            
            if(DN)
            {
                DN->updateGeometry(simulationParameters.dt);
                updateLoadControllers(simulationParameters.runID, false);
                
                DN->assembleAndSolveGlide(simulationParameters.runID);
                simulationParameters.dt=DDtimeIntegrator<0>::getGlideTimeIncrement(*DN); // TO DO: MAKE THIS std::min between DN and CrackSystem

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////CAMERON'S TIME EDIT///////////////////////////////////////////////////////////////////////////////////////
		simulationParameters.dt=1*pow(10,0);
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
                // output
                DN->io().output(simulationParameters.runID);

                
                //                for(const auto& loop : DN->loops())
                //                {
                //                    if(loop.second->loopType==DislocationLoopIO<dim>::GLISSILELOOP)
                //                    {
                //                        PlanarDislocationSuperLoop<typename DislocationNetworkType::LoopType> superLoop(*loop.second);
                //                    }
                //                }
                
                // move
                DN->moveGlide(simulationParameters.dt);
                
                // menage discrete topological events
                DN->singleGlideStepDiscreteEvents(simulationParameters.runID);
            }
            simulationParameters.totalTime+=simulationParameters.dt;
            ++simulationParameters.runID;
        }
        
        /**********************************************************************/
        void runGlideSteps()
        {/*! Runs a number of simulation time steps defined by simulationParameters.Nsteps
          */
            const auto t0= std::chrono::system_clock::now();
            while (simulationParameters.runID<simulationParameters.Nsteps)
            {
                model::cout<<std::endl; // leave a blank line
                singleGlideStep();
            }
            model::cout<<greenBoldColor<<std::setprecision(3)<<std::scientific<<simulationParameters.Nsteps<< " simulation steps completed in "<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" [sec]"<<defaultColor<<std::endl;
        }
        
        /**********************************************************************/
#ifdef _MODEL_GREATWHITE_
#include <DefectiveCrystalGreatWhite.h>
#endif
        
        
        /**********************************************************************/
        VectorDim displacement(const VectorDim& x) const
        {/*!\param[in] P position vector
          * \returns The displacement field in the DefectiveCrystal at P
          */
            VectorDim temp(VectorDim::Zero());
            if(DN)
            {
                temp+=DN->displacement(x);
            }
            if(CS)
            {
                temp+=CS->displacement(x);
            }
            return temp;
        }
        
        /**********************************************************************/
        void displacement(std::vector<FEMnodeEvaluation<ElementType,dim,1>>& fieldPoints) const
        {
            if(DN)
            {
                DN->displacement(fieldPoints);
            }
            if(CS)
            {
                CS->displacement(fieldPoints);
            }
        }
        
        /**********************************************************************/
        MatrixDim stress(const VectorDim& x) const
        {/*!\param[in] P position vector
          * \returns The stress field in the DefectiveCrystal at P
          * Note:
          */
            MatrixDim temp(MatrixDim::Zero());
            if(DN)
            {
                temp+=DN->stress(x);
            }
            if(CS)
            {
                temp+=CS->stress(x);
            }
            return temp;
        }
        
        /**********************************************************************/
        MatrixDim plasticDistortion() const
        {/*!\param[in] P position vector
          * \returns The stress field in the DefectiveCrystal at P
          * Note:
          */
            MatrixDim temp(MatrixDim::Zero());
            if(DN)
            {
                temp+=DN->plasticDistortion();
            }
            if(CS)
            {
                temp+=CS->plasticDistortion();
            }
            return temp;
        }
        
        /**********************************************************************/
        MatrixDim plasticDistortionRate() const
        {/*!\param[in] P position vector
          * \returns The stress field in the DefectiveCrystal at P
          * Note:
          */
            MatrixDim temp(MatrixDim::Zero());
            if(DN)
            {
                temp+=DN->plasticDistortionRate();
            }
            if(CS)
            {
                temp+=CS->plasticDistortionRate();
            }
            return temp;
        }
        
        /**********************************************************************/
        MatrixDim plasticStrainRate() const
        {/*!\param[in] P position vector
          * \returns The stress field in the DefectiveCrystal at P
          * Note:
          */
            MatrixDim temp(plasticDistortionRate());
            return 0.5*(temp+temp.transpose());
        }
        
    };
}
#endif
