#include <iostream>
#include <cmath>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <iomanip>
#include <time.h>
#include <vector>

#include <iomanip>
#include <iostream>
#include <fstream>
#include<map>

#include <Vacancies_functions.h>

#ifndef Test_H_
#define Test_H_

using namespace std;


class DislocationStructure{
//Class to hold the dislocation node positions, loops, and node connectivities. Also included are nodal velocities (typically 0 for this application), burgers vector, plane normal, and boundary information
public:
	//No constructor needed at the moment

	std::map<int, std::array<double, 6>> nodePositions; //NODES Key: nodeID, value: node position [0-2], nodal velocities [3-5]
	std::map<int, std::array<double, 6>> loops; //LOOPS Key: loopID, value[0-2]: burgers vector ,value[3-6]: plane normal
	std::map<std::array<int, 2>, std::array<int, 2>> nodeConnectivity; //CONNECTIVITY Key: Connectivity, value[0] = loopnumber, value[1] = bool for if it is a boundary link or not (1=boundary link, 0=internal link

	int meshType; //Defines what geometry of mesh we are using. meshType=0 is a cylindrical mesh, meshType=1 is a prisim mesh
	double height, radius;
	double xWidth, yWidth; 
	//For meshType=0 (cylinder)

	//For meshType=1 (rect prisim) - radius stays undefined


	void readNodesAndSegments(int runID); //Function to update the nodes and segments 
	void createNodesAndSegments(int numNodes, double z_dist, double x_dist); //Function to create the nodes and segments from pre-defined instructions
	void printevl(int runID); //Function to print the required evl file for VTK visualization

};

//void readNodesAndSegments(int runID,std::map<int,std::array<double,3>>& nodePositions, std::map<int,std::array<double, 6>>& loops, std::map<std::array<int, 2>,std::array<int, 2>>& nodeConnectivity)
void DislocationStructure::readNodesAndSegments(int runID)
//Function that generates an array with all of the nodesNums - their current positions
{

	std::cout << "Updating DislocationStructure..." << std::endl;

	int IDnum = runID-1;
	std::string fileString= "./evl/evl_" + std::to_string(IDnum) + ".txt";  //Creates file name based on the current runID
	//std::cout << fileString << std::endl;

	ifstream ifEVLStream;
	ifEVLStream.open(fileString);

	ifEVLStream.precision(15);

	int numNodes, numLoops, numLoopSegments;

	ifEVLStream >> numNodes;
	ifEVLStream >> numLoops;
	ifEVLStream >> numLoopSegments;

	//std::cout << numNodes << std::endl;
	//std::cout << numLoops << std::endl;
	//std::cout << numLoopSegments << std::endl;

	///////////////////////
	//NODES////////////////
	///////////////////////
	int nodeID;
	std::array<double,6> position;

	for(int i=0; i<numNodes;i++)
		{
		ifEVLStream >> nodeID;
		ifEVLStream >> std::setprecision(15) >>position[0];//All the nodal positions
		ifEVLStream >> std::setprecision(15) >>position[1];
		ifEVLStream >> std::setprecision(15) >>position[2];

		ifEVLStream >> std::setprecision(15) >>position[3];//All the nodal velocities
		ifEVLStream >> std::setprecision(15) >>position[4];
		ifEVLStream >> std::setprecision(15) >>position[5];

		//std::cout << nodeID << "     "; 
		//std::cout << std::setprecision(15) << std::scientific << position[0] << "  "; 
		//std::cout << position[1] << "  "; 
		//std::cout << position[2] << std::endl;
		
		nodePositions.insert( std::make_pair(nodeID, position)  );
		//nodePositions[nodeID] = position;
		ifEVLStream.ignore(256, '\n');
		}

	///////////////////////
	//LOOPS////////////////
	///////////////////////
	int loopID;
	std::array<double,6> burgersAndTangent;

	for(int i=0; i<numLoops;i++)
		{
		ifEVLStream >> loopID;
		for(int j=0; j<6;j++)
			ifEVLStream >> burgersAndTangent[j];
		
		//std::cout << loopID << "     "; 

		//for(int j=0; j<6;j++)
			//std::cout<< burgersAndTangent[j] << "   ";

		//std::cout << std::endl;

		loops.insert( std::make_pair(loopID, burgersAndTangent)  );
		ifEVLStream.ignore(256, '\n');
		}

	///////////////////////////
	//LOOP SEGMENT CONNECTIVITY
	///////////////////////////
	std::array<int,2> connectivity;
	std::array<int,2> nodeInfo;

	for(int i=0; i<numLoopSegments;i++)
		{
		ifEVLStream >> nodeInfo[0]; //loop ID of segment
		
		ifEVLStream >> connectivity[0]; //node 1 of segment
		ifEVLStream >> connectivity[1]; //node 2 of segment

		ifEVLStream >> nodeInfo[1]; //whether or not it lies on a boundary

		//std::cout << nodeInfo[0] << "   " << connectivity[0] << "   " << connectivity[1] << "   " << nodeInfo[1] << std::endl;

		nodeConnectivity.insert( std::make_pair(connectivity, nodeInfo)  );
		ifEVLStream.ignore(256, '\n');
		}

	ifEVLStream.close();

}

void DislocationStructure::createNodesAndSegments(int numNodes, double z_dist, double x_dist)
//Function to create a dislocation structure from a given number of nodes (distributed across the straight dislocation) and an x and y overhang
{

	///////////////////////
	//NODES////////////////
	///////////////////////
	int nodeID;
	std::array<double,6> position;

	int totalNodes = numNodes+2;

	double totalHeight = L3*(1+z_dist);
	double totalWidth = L1*(1+x_dist);

	for(int i=0; i<totalNodes;i++)
		{
		nodeID = i;

		if(i<numNodes)
		{
			position[0]=0;//All the nodal positions
			position[1]=0;
			position[2]=-totalHeight/2+totalHeight*i/(numNodes-1);
		}

		if(i==numNodes)
		{
			position[0]=-totalWidth;//All the nodal positions
			position[1]=0;
			position[2]=totalHeight/2;
		}

		if(i==numNodes+1)
		{
			position[0]=-totalWidth;//All the nodal positions
			position[1]=0;
			position[2]=-totalHeight/2;
		}

		position[3]=0;//All the nodal velocities
		position[4]=0;
		position[5]=0;

		
		nodePositions.insert( std::make_pair(nodeID, position)  );

		}

	///////////////////////
	//LOOPS////////////////
	/////////////////////// ....1 loop
	int loopID;
	std::array<double,6> burgersAndTangent;

	loopID=0;

	burgersAndTangent[0]=0; //Burgers vector
	burgersAndTangent[1]=1;
	burgersAndTangent[2]=0;

	burgersAndTangent[3]=1; //Plane normal
	burgersAndTangent[4]=0;
	burgersAndTangent[5]=0;

	loops.insert( std::make_pair(loopID, burgersAndTangent)  );

	///////////////////////////
	//LOOP SEGMENT CONNECTIVITY
	///////////////////////////
	std::array<int,2> connectivity;
	std::array<int,2> nodeInfo;

	int numLoopSegments = 3+numNodes-1;

	for(int i=0; i<numLoopSegments;i++)
		{
		nodeInfo[0]=i; //loop ID of segment
		
		connectivity[0]=i; //node 1 of segment
		connectivity[1]=i+1; //node 2 of segment
			

		if(i==numLoopSegments-1)
			connectivity[1]=0; //Making sure the last segment connects from the final node to the first node

		if(i<numNodes-1)		
			nodeInfo[1]=0; //whether or not it lies on a boundary (0=not on boundary, 1 = on boundary) -> determines whether segment will be VISIBLE or not
		else
			nodeInfo[1]=1; //boundary node!
		

		//std::cout << nodeInfo[0] << "   " << connectivity[0] << "   " << connectivity[1] << "   " << nodeInfo[1] << std::endl;

		nodeConnectivity.insert( std::make_pair(connectivity, nodeInfo)  );
		}


}


void DislocationStructure::printevl(int IDnum)
//Function that prints out an evl file
{
	
	if(outputEVL==1)	
	{

		std::cout << "Printing evl_" << IDnum << std::endl;

		//int IDnum = runID;
		std::string fileString= "./evl1/evl_" + std::to_string(IDnum) + ".txt";  //Creates file name based on the current runID
		
		ofstream OutEVLStream;
		OutEVLStream.open(fileString);

		//First pint number of nodes, number of loops, and number of loop segments

		OutEVLStream << nodePositions.size() << "\n";	//# nodes
		OutEVLStream << loops.size() << "\n";		//# loops
		OutEVLStream << nodeConnectivity.size() << "\n";//# loop segments

		double ZeroVal=0;
		double OneVal =1;

		//Print all of the node positions, velocities, and other relevant data
		for( auto const& [key, val] : nodePositions )
		{
			OutEVLStream << key << "\t";
			
			OutEVLStream << std::scientific<< std::setprecision(15) << std::right << nodePositions[key][0] << " " << nodePositions[key][1] << " " << nodePositions[key][2] << "\t";
			OutEVLStream << std::scientific << std::setprecision(15) << std::right << ZeroVal << " " << ZeroVal << " " << ZeroVal << "\t" << OneVal << "\t";
			OutEVLStream << 0 << "\t" << 0 << "\t" << key << "\n";
		}

		//Print all of the node loops
		for( auto const& [key, val] : loops )
		{
			OutEVLStream << key << "\t";
			
			OutEVLStream << std::scientific << std::setprecision(15) << std::right << loops[key][0] << " " << loops[key][1] << " " << loops[key][2] << "\t";
			OutEVLStream << std::scientific << std::setprecision(15) << std::right << loops[key][3] << " " << loops[key][4] << " " << loops[key][5] << "\t";

			OutEVLStream << std::scientific << std::setprecision(15) << std::right << ZeroVal << " " << ZeroVal << " " << ZeroVal << "\t";
			OutEVLStream << 0 << "\t" << 0 << "\t";
			OutEVLStream << ZeroVal << "\t" << ZeroVal << "\t" << ZeroVal <<"\n";
		}

		//Print all of the node connectivities
		for( auto const& [key, val] : nodeConnectivity)
		{
			OutEVLStream << 0 << "\t" << key[0] << "\t" << key[1] << "\t" << val[1] << "\n"; // <- This accounts for segments being boundary or non-boundary segments
			//OutEVLStream << 0 << "\t" << key[0] << "\t" << key[1] << "\t" << 0 << "\n"; // <- This does NOTE account for segments being boundary or non-boundary segments
		}

		OutEVLStream.close();


		//Print the reqired ddAux_X.txt file needed for DDvtk
		std::cout << "Printing ddAux_" << IDnum << std::endl;

		//int IDnum = runID;
		std::string fileStringddAux= "./evl1/ddAux_" + std::to_string(IDnum) + ".txt";  //Creates file name based on the current runID
		
		ofstream OutddAuxStream;
		OutddAuxStream.open(fileStringddAux);

		OutddAuxStream << 0 << "\n" << 0 << "\n" << 0;

		OutddAuxStream.close();

	}

}


double returnDistanceToSegment(double segmentNode1[3], double segmentNode2[3], double VacancyPoint[3])
//Function to return the shortest distance between the vacancy point and the line segment
//Can follow from http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
{

	double temp1[3], temp2[3], temp3[3]; 
	double temp1Mag, temp2Mag, temp3Mag; //Magnitude scalars
	double numeratorCrossProduct[3], SegmentVec[3];

	double theta1, theta2;

	double numeratorMag, denominatorMag;

	double d; //final computed distance from dislocation segment

	//std::cout << std::endl;
	//std::cout << "Node 1:" << segmentNode1[0] << "  " << segmentNode1[1] << "  " << segmentNode1[2] << std::endl;
	//std::cout << "Node 2:" << segmentNode2[0] << "  " << segmentNode2[1] << "  " << segmentNode2[2] << std::endl;
	//std::cout << "Vacancy:" << VacancyPoint[0] << "  " << VacancyPoint[1] << "  " << VacancyPoint[2] << std::endl;

	for(int i = 0;i<3;i++)
	{
		temp1[i] = VacancyPoint[i]-segmentNode1[i];
		temp2[i] = VacancyPoint[i]-segmentNode2[i];
		temp3[i] = segmentNode1[i]-segmentNode2[i];
	}

	//std::cout << "Temp1 = " << temp1[0] << " " << temp1[1] << " " << temp1[2] << std::endl;
	//std::cout << "Temp2 = " << temp2[0] << " " << temp2[1] << " " << temp2[2] << std::endl;
	//std::cout << "Temp3 = " << temp3[0] << " " << temp3[1] << " " << temp3[2] << std::endl;
	//Compute the angles between the vectors to see if the point lays in or out of the segment

	temp1Mag  = pow(pow(temp1[0],2) + pow(temp1[1],2) + pow(temp1[2],2),0.5);
	temp2Mag  = pow(pow(temp2[0],2) + pow(temp2[1],2) + pow(temp2[2],2),0.5);
	temp3Mag  = pow(pow(temp3[0],2) + pow(temp3[1],2) + pow(temp3[2],2),0.5); //Calculate magnitudes 

	//std::cout << "Temp1mag = " << temp1Mag << std::endl;
	//std::cout << "Temp2mag  = " << temp2Mag << std::endl;
	//std::cout << "Temp3mag  = " << temp3Mag << std::endl;


	//Using law of cosines to find the angle between segments - used to determine projection of vacancy point on segment line
	theta1 = acos( ( pow(temp3Mag,2) + pow(temp2Mag,2) - pow(temp1Mag,2) )/(2*temp3Mag*temp2Mag) )*180/3.14159;
	theta2 = acos( ( pow(temp3Mag,2) - pow(temp2Mag,2) + pow(temp1Mag,2) )/(2*temp3Mag*temp1Mag) )*180/3.14159;

	//std::cout << "Theta 1 = " << theta1 << std::endl;
	//std::cout << "Theta 2 = " << theta2 << std::endl;

	//Only project the point onto the line if it lays within th 3D space projected orthogonally to the line segment
	if ((theta1<90) && (theta2<90)) 
	{
		//std::cout << "Method 1" << std::endl;

		//Compute the cross product of the numerators (temp1 and temp2)
		numeratorCrossProduct[0] = (temp1[1]*temp2[2] - temp1[2]*temp2[1]);
		numeratorCrossProduct[1] = -(temp1[0]*temp2[2] - temp1[2]*temp2[0]);
		numeratorCrossProduct[2] = (temp1[0]*temp2[1] - temp1[1]*temp2[0]);

		//Find magnitudes
		numeratorMag = pow(pow(numeratorCrossProduct[0],2) + pow(numeratorCrossProduct[1],2) + pow(numeratorCrossProduct[2],2),0.5);
		denominatorMag = pow(pow(temp3[0],2) + pow(temp3[1],2) + pow(temp3[2],2),0.5);
		
		
		d = abs(numeratorMag/denominatorMag);
	}

	else //If outside of the linear projection of the line, find the closest node on the segment
	{
		//std::cout << "Method 2" << std::endl;
	 	d = min(pow(pow(temp1[0],2) + pow(temp1[1],2) + pow(temp1[2],2),0.5), pow(pow(temp2[0],2) + pow(temp2[1],2) + pow(temp2[2],2),0.5));
	}

	//std::cout << "Distance is " << d << std::endl;

	return d;
}

void findClosestSegment(double position[3], DislocationStructure DS, double distanceInfo[])
//Function to find the distance to the closest segment along with the 2 node IDs needed to identify the segment
//the array distanceInfo is an array with size 3 where:
// distanceInfo[0] = distance to segment in b units
// distanceInfo[1] = nodeID1, 
// distanceInfo[2] = nodeID2 , 
//note that the nodeID's are left in double form, but should be converted to integers
{

	distanceInfo[0] =100000000; //Distance of closest segment set to an arbitrarily high value
	double tempd;
	
	double tempPos1[3] ,tempPos2[3], positionInB[3];

	//Convert our absolute position system to the burger's position system

	for(int i = 0;i<3;i++)
		positionInB[i] = position[i]/b;
		

	for( auto const& [key, val] : DS.nodeConnectivity )
	{
		if (val[1]==0) //Ensures that we are looking at a non-boundary segment
			{
			for(int i =0;i<3;i++)
				{
					tempPos1[i] = DS.nodePositions[key[0]][i];
					tempPos2[i] = DS.nodePositions[key[1]][i];	
				}
		
			tempd = returnDistanceToSegment(tempPos1, tempPos2, positionInB);
			
			if(tempd<distanceInfo[0])
				{
					distanceInfo[0] = tempd;
					distanceInfo[1] = key[0];
					distanceInfo[2] = key[1];
				}
			}

		//std::cout << "key = " << key[0] << "  " << key[1] << std::endl;

	}

	//std::cout << "Shortest distance is " << distanceInfo[0] << std::endl;
	//std::cout << "Node " << distanceInfo[1] << " to node " << distanceInfo[2] << std::endl;

}

void printStatistics(DislocationStructure DisStruct, long int runID)
//Function to print the running values of vacancies absorbed, dislocation velocity, timestep, and other relevant statistics
{

	ofstream OutStatsStream;
	OutStatsStream.open("./results/stats.txt", std::ios_base::app); //Stats file to hold the calculated data from each output -- the stats are calculated every PrintStatsFreq trial, defined in Vacancy_Inputs

	int incrementalVacsAbsorbed = RunningVacAbsorbed - RunningLastVacNumber; //Calculate the number of vacancies absorbed since the previous calculation
	RunningLastVacNumber = RunningVacAbsorbed; //Update the counter for the number of vacancies absorbed

	int incrementalVacsEmitted = RunningVacEmitted - RunningLastVacEmitted; //Calculate the number of vacancies emitted since the previous calculation
	RunningLastVacEmitted = RunningVacEmitted; //Update the counter for the number of vacancies emitted


	double incrementalTimeStep = totalGlobalTime - lastTotalGlobalTime; //How much time has passed since last calculation
	lastTotalGlobalTime = totalGlobalTime; //Reset the marker for the last global time calculated


	/*
	double totalAreaMoved = 0;

	//Loop through the segments to extract the area swept out -> calculate distance traveled
	for( auto const& [key, val] : DisStruct.nodeConnectivity )
	{
		if(val[1]==1) //Ignore the segment if it is a boundary segment
			continue;

		double x_dif = DisStruct.nodePositions[key[0]][0] - DisStruct.nodePositions[key[1]][0];  
		double z_dif = abs(DisStruct.nodePositions[key[0]][2] - DisStruct.nodePositions[key[1]][2]);

		if(x_dif>0) //Ignore all screw connecting segments
			continue;

		//totalAreaMoved = totalAreaMoved + DisStruct.nodePositions[key[0]][0]*z_dif; //Calculate the total area swept out the segment
		totalAreaMoved = totalAreaMoved + DisStruct.nodePositions[key[0]][0]*z_dif; //Calculate the total area swept out the segment
	}
	*/

	////////Calculate the total and incremental distance moved by using the number of vacancies absorbed and emitted//////
	double segmentLength = L1/(numNodes-1); //Individual dislocation segment length in b units

	double atomicVolumeinb = atomicVolume*pow(b,-3); //Atomic volume in b units

	double h = ( atomicVolumeinb/ segmentLength); //Climb height in b units according to volume swept out by absorption of one vacancy in one segemnet

	double totalDistanceMoved = (RunningVacAbsorbed-RunningVacEmitted)*h/(numNodes-1); //Normalize the total area swept out to the entire length of the segment

	//double totalDistanceMoved = totalAreaMoved/L1; //Normalize the total area swept out to the entire length of the segment

	double totalDislocationVelocity = totalDistanceMoved*b/totalGlobalTime; //Calculate the total average dislocation velocity in ABSOLUTE UNITS

	//Calculate instantaneous dislocation velocity
	double incrementalDislocationDistance = (incrementalVacsAbsorbed - incrementalVacsEmitted)*h/(numNodes-1);;	
	double incrementalDislocationVelocity = incrementalDislocationDistance*b/incrementalTimeStep; //Incremental dislocation velocity in ABSOLUTE UNITS

	//lastDistanceMoved = totalDistanceMoved; //Update the last distance moved


	////////Print all statistics consistent across ALL trials////////

	//Print all the variables describing the simulation state (temp, concentration, runID, applied stress, etc.)//
	OutStatsStream << Temp << " " << vacancyConcentration << " " << runID << " ";

	if(useStress==1) //If an applied stress is used, print the stress values
		OutStatsStream << appliedPressure << " ";

	if(dynamicBoxResizing==1)//Output the L1, L2, L3 sizes to see side lengths if dynamicresizing is on
		OutStatsStream << L1 << " " << L2 << " " << L3 << " ";
		

	//Print all the total # vacancies absorbed, # vacancies absorbed during the last step, total global time, time since last printing, total dislcoation velocity, and incremental dislocation velocity
	OutStatsStream << RunningVacAbsorbed << " " << incrementalVacsAbsorbed << " "; 
	OutStatsStream << RunningVacEmitted << " " << incrementalVacsEmitted << " "; 
	OutStatsStream << totalGlobalTime << " " << incrementalTimeStep << " " << totalDislocationVelocity << " " << incrementalDislocationVelocity << " ";
	
	OutStatsStream << std::endl;
	
	OutStatsStream.close();


	//////Print a header file at the beginning of each trial//////
	if( (Temp==minTemp) && (RunningVacAbsorbed == incrementalVacsAbsorbed) && (runID<(PrintStatsFreq+2)) ) //Reprint the header file during the initial simulations
	{
		ofstream OutStatsHeaderStream;
		OutStatsHeaderStream.open("./results/stats_header.txt"); //File to hold the header labels for the stats writeup

		OutStatsHeaderStream << "Temperature [K]" << std::endl;
		OutStatsHeaderStream << "Vacancy_Concentration [vacs/m^3]" << std::endl;
		OutStatsHeaderStream << "runID" << std::endl;

		if(useStress==1) //If an applied stress is used
			OutStatsHeaderStream << "Pressure [Pa]" << std::endl;

		if(dynamicBoxResizing==1)//Output the L1, L2, L3 sizes to see side lengths if dynamicresizing is on
			OutStatsHeaderStream << "L1 [b]" << std::endl << "L2 [b]" << std::endl << "L3 [b]" << std::endl;

		OutStatsHeaderStream << "Total_Vacancies_Absorbed [#_vac]" << std::endl;
		OutStatsHeaderStream << "Incremental_Vacancies_Absorbed [#_vac]" << std::endl;
		OutStatsHeaderStream << "Total_Vacancies_Emitted [#_vac]" << std::endl;
		OutStatsHeaderStream << "Incremental_Vacancies_Emitted [#_vac]" << std::endl;
		OutStatsHeaderStream << "Total_Time [s]" << std::endl;
		OutStatsHeaderStream << "Incremental_Time [s]" << std::endl;
		OutStatsHeaderStream << "Total_Avg_Dislocation_Velocity [m/s]" << std::endl;
		OutStatsHeaderStream << "Incremental_Dislocation_Velocity [m/s]" << std::endl;

		OutStatsHeaderStream.close();
	}
	///////////////////////////////
}


void remeshEmission(DislocationStructure& DisStruct)
//Function to remesh the dislocation structure after an emission event
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

			std::cout << "Node " << key[0] << " and " << key[1] << std::endl;

			double x_dif = DisStruct.nodePositions[key[0]][0] - DisStruct.nodePositions[key[1]][0];  
			double y_dif = DisStruct.nodePositions[key[0]][1] - DisStruct.nodePositions[key[1]][1];
			double z_dif = DisStruct.nodePositions[key[0]][2] - DisStruct.nodePositions[key[1]][2];

			x_dif = abs(x_dif);
			y_dif = abs(y_dif);
			z_dif = abs(z_dif);


			int firstNode = key[0];
			int secondNode = key[1];

			std::cout << "x_dif = " << x_dif << std::endl;
			std::cout << "y_dif = " << y_dif << std::endl;
			std::cout << "z_dif = " << z_dif << std::endl;


			//First check to see if there is an angled segment, if so, add a node to make it straight
			if (x_dif>0 && z_dif>0)
			{

				if(key[0]==key[1])
					{
					continue; //Sanity Check
					}

				std::cout << "Looping at curvature node" << std::endl;

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

				std::cout << "(1)New node num = " << newNodeNum << std::endl;
				std::cout << "(1)key[1] = " << secondNode << std::endl;

				//Create and insert the new node
				std::array<double,6> position;

				
				//Select the HIGHEST x-value of the two
				position[0]= max(DisStruct.nodePositions[key[0]][0], DisStruct.nodePositions[key[1]][0]); //Max x-value 
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

				std::cout << "Curvature between node "<< key[0] << " at :" << DisStruct.nodePositions[key[0]][0] << ", " << DisStruct.nodePositions[key[0]][1] << ", " << DisStruct.nodePositions[key[0]][2] << std::endl;
				std::cout << " and node " << key[1] << " at " << DisStruct.nodePositions[key[1]][0] << ", " << DisStruct.nodePositions[key[1]][1] << ", " << DisStruct.nodePositions[key[1]][2] << std::endl;
				std::cout << "Adding node " << newNodeNum << " at : " << position[0] << ", " << position[1] << ", " << position[2] << std::endl;

				std::cout << "(2)New node num = " << newNodeNum << std::endl;
				std::cout << "(2)key[1] = " << secondNode << std::endl;

				std::cout << "Inserting node " << newNodeNum << std::endl;
				DisStruct.nodePositions.insert( std::make_pair(newNodeNum, position)  );

				//Update the node connectivity 

				//Delete old connectivity
				std::array<int, 2> nodes;

				nodes[0] = key[0];
				nodes[1] = key[1];

				std::cout << "Deleting segment " << nodes[0] <<"->" <<nodes[1] <<std::endl;
				DisStruct.nodeConnectivity.erase(nodes); //erase the old connectivity value

				//Insert the new connectivity values
				std::array<int,2> connectivity;
				//connectivity[0] = key[0];
				connectivity[0] = firstNode;
				connectivity[1] = newNodeNum;

				std::cout << "(3)New node num = " << newNodeNum << std::endl;
				std::cout << "(3)key[1] = " << secondNode << std::endl;

				assert(key[0]!=newNodeNum); //Sanity Check

				std::array<int,2> nodeInfo;
				nodeInfo[0] = 0;
				nodeInfo[1] = 0;

				std::cout << "Adding segment " << connectivity[0] <<"->" <<connectivity[1] <<std::endl;
				DisStruct.nodeConnectivity.insert( std::make_pair(connectivity, nodeInfo)  );

				//Second connectivity value
				std::array<int,2> connectivity1;
				connectivity1[0] = newNodeNum;
				//connectivity[1] = key[1];
				connectivity1[1] = secondNode;

				std::cout << "(4)New node num = " << newNodeNum << std::endl;
				std::cout << "(4)key[1] = " << secondNode << std::endl;

				assert(secondNode!=newNodeNum); //Sanity Check

				std::cout << "Adding segment " << connectivity1[0] <<"->" <<connectivity1[1] <<std::endl;
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

				std::cout << "Overlapping nodes found" << std::endl;

				std::cout << "Node " << key[0] << " :";
				for(int i = 0; i <3; i++)
					{
					std::cout << DisStruct.nodePositions[key[0]][i];
					}
				
				std::cout << std::endl;

				std::cout << "Node " << key[1] << " :";
				for(int i = 0; i <3; i++)
					{
					std::cout << DisStruct.nodePositions[key[1]][i];
					}
										
				std::cout << std::endl;

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
				//connectivity[0] = key[0];
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

template <typename Vacancy>
void emissionEvent(vector<Vacancy> &vacancyArray, DislocationStructure& DisStruct)
//Function to facilitate the negative climb and addition of a new vacancy to the simulation box
{

	int node1, node2;
	////Select a random dislocation segment and have it climb in the negative direction - only select edge segments!////

	double probOfSelection = pow(numNodes-1, -1); //Probability of one of the straight edge segements being selected to emit a vacancy
	
	int isEmitted = 0; //Int for signifying that emission has occured

	while(isEmitted==0) //Continue to loop through until emission has occured
		{

		for( auto const& [key, val] : DisStruct.nodeConnectivity) //Loop through the dislocation segments and randomly select an edge segment
				{

					if(key[0]==0 || key[1]==(numNodes-1) ) //Skip emission from the first or last segment, for now. 
						{
						std::cout << " Test 1" << std::endl;
						continue;
						}

					if(key[0]==key[1])
						{
						std::cout << " Test 2" << std::endl;
						continue; //Sanity Check
						}

					double z_dif = abs(DisStruct.nodePositions[key[0]][2] - DisStruct.nodePositions[key[1]][2]);

					if( (z_dif==0) || (val[1]==1) ) //Skip through until a proper edge segment is randomly selected
						{
						std::cout << " Test 3" << std::endl;
						continue;
						}

					double randomChance = ZeroOnedistribution(generator);
					std::cout << "Selection prob = " << probOfSelection << " and chance of " << randomChance << std::endl;

					if( probOfSelection>randomChance ) //Select a segment with a uniform probability across all straight edge segments
						{

						node1 = key[0];
						node2 = key[1];

						std::cout << "Selected node " << key[0] << " and node " << key[1] << " to negatively climb" << std::endl;


						std::cout << "Node " << key[0] << " :";
						for(int i = 0; i <3; i++)
							{
							std::cout << DisStruct.nodePositions[key[0]][i] << " ";
							}
						
						std::cout << std::endl;

						std::cout << "Node " << key[1] << " :";
						for(int i = 0; i <3; i++)
							{
							std::cout << DisStruct.nodePositions[key[1]][i] << " ";
							}

						std::cout << std::endl;

						isEmitted = 1; //Set isEmitted to 1 so the while loop will be exitted

						break; //Exit the loop once a segment has been selected

						}
					else
						{
						std:cout <<"Nodes skipped " << std::endl;
						}

				}
		}	

		
		/////////////////////////////////////////////////////////////////////////////////////////				
		//Execute NEGATIVE CLIMB of the segment that the vacancy was absorbed into - via movement of the two segment
		/////////////////////////////////////////////////////////////////////////////////////////

		//double atomicVolume = (4/3)*3.14159*pow(atomicRadius/b, 3); //Atomic volume in b units
		//Use the pre-determined atomic volume
		double segmentLength = L1/(numNodes-1); //Individual dislocation segment length in b units

		double atomicVolumeinb = atomicVolume*pow(b,-3); //Atomic volume in b units

		double h = ( atomicVolumeinb/ segmentLength); //Climb height in b units according to volume swept out by absorption of one vacancy

		//DisStruct.nodePositions[node1][0]+=-h;

		//DisStruct.nodePositions[node2][0]+=-h;

		DisStruct.nodePositions[node1][0]+=-10;

		DisStruct.nodePositions[node2][0]+=-10;
		/////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////
		

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
	////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////
	remeshEmission(DisStruct); //Function to remesh the dislocation structure after an emission event
	////////////////////////////////////////////////////////////////////
}


template <typename Vacancy>
void vacancyEmission(vector<Vacancy> &vacancyArray, DislocationStructure& DisStruct, double timeStep)
//Function to randomly emit vacancies from the dislocation segments based on temperature, stress, and the global vacancy concentration
{
	double dislocationDensity = 1/(L1*L2*pow(b,2)); //Dislocation density [m^-2]
	double simVolume = L1*L2*L3*pow(b,3); //Total mesh volume
	int numVacancies = static_cast<int>(vacancyArray.size()); //Number of vacancies currently in the mesh
	double currentVacancyConcentration = numVacancies/simVolume; //Current calculated vacancy concentration

	double formationEnthalpy = vacFormationEnergy - appliedPressure*(1+volumetricStrain); //Stress-dependent vacancy formation energy
	
	double thermalVacConcentration = (density*Na*1*exp(-1*formationEnthalpy/(Temp*Kb))/molarMass ) *pow(100,3); //Thermal concentration of vacancies considering current stress state

	double emissionRate = 2*3.1415* simVolume * dislocationDensity * v*exp(-Eo/(Kb*Temp)) * (1 - currentVacancyConcentration/thermalVacConcentration) / b; //Calculation of the emission rate!

	std::cout << "Checking for Vacancy Emissions..." << std::endl;
	std::cout << "Thermal Vacancy Concentration [ideal] = " << thermalVacConcentration << " vacs/m^3" << std::endl;
	std::cout << "Actual Vacancy Concentration = " << currentVacancyConcentration << " vacs/m^3" << std::endl;
	std::cout << "Emission rate = " << emissionRate << " vacs/s" << std::endl;

	if (emissionRate>0)
		{
			std::cout << "Timestep = " << timeStep << " --> Emission probability = " << emissionRate*timeStep*100 << "%" << std::endl;

			if (emissionRate*timeStep>ZeroOnedistribution(generator)) //Only consider vacancy emission if the event is randomly selected
			{
				RunningVacEmitted++;
				std::cout <<"1 vacancy emitted" << std::endl;

				emissionEvent(vacancyArray,DisStruct);  //Perform the negative climb and new vacancy addition
				
			}
			else
				std::cout <<"No vacancy emitted" << std::endl;
	
		}
	else
 		std::cout << "Timestep = " << timeStep << " --> Emission probability = 0" << "%" << std::endl;


	std::cout <<"Total Vacancies Emitted = " << RunningVacEmitted << std::endl;

}


#endif

