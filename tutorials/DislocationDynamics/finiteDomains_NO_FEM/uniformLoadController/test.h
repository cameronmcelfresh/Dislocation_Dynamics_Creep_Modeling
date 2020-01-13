#include <iostream> 
#include <string>
#include <cstdio>
#include <fstream>

#ifndef HEADER1_H
#define HEADER1_H


void rewriteTemperature(double newTemp)
//Function to rewrite the W material input file with the updated temperature value -- brute force method used
//Temperature is given in Kelvin - current function is designed to update the W material input file
{

	using namespace std;

	//Create a text file to hold the new W info called W_new, the delete the old W file and change the name of the new file

	/*
	ofstream OutMaterialStream;
	OutMaterialStream.open("../../MaterialsLibrary/W.txt"); //W_new file to hold new W info

	OutMaterialStream << "materialName=W;\n";
	OutMaterialStream << "crystalStructure=BCC;\n";
	OutMaterialStream << "b_SI=0.2722e-9;	# [m] 		Burgers vector magnitude\n";
	OutMaterialStream << "mu0_SI=161.0e9;	# [Pa] 		temperature-independent shear modulus coeff in mu=mu0+mu1*T\n";
	OutMaterialStream << "mu1_SI=0.0;		# [Pa/K] 	temperature-dependent shear modulus coeff in mu=mu0+mu1*T\n";
	OutMaterialStream << "nu=0.28;		# [-]		Poisson's ratio\n";
	OutMaterialStream << "rho_SI=19250.0;	# [kg/m^3]	mass density\n";
	OutMaterialStream << "T=" <<newTemp <<";			# [K]		temperature\n";
	OutMaterialStream << "Tm=3695.0;		# [K]		melting temperature\n";
	OutMaterialStream << "\n";

	OutMaterialStream << "# Mobility parameters (G.Po et al. A phenomenological dislocation mobility law for bcc metals. Acta Mater 119, 2016)\n";
	OutMaterialStream << "B0e_SI=4.26e-04;	# [Pa*s]	temperature-independent drag coefficient for edge  dislocations\n";
	OutMaterialStream << "B1e_SI=0.87e-06;	# [Pa*s/K]	temperature-dependent drag coefficient for edge  dislocations\n";
	OutMaterialStream << "B0s_SI=9.8e-4;		# [Pa*s]	temperature-independent drag coefficient for screw  dislocations\n";
	OutMaterialStream << "B1s_SI=0.0;			# [Pa*s/K]	temperature-dependent drag coefficient for screw  dislocations\n";
	OutMaterialStream << "Bk_SI=8.3e-05;		# [Pa*s]	drag coefficient for kinks\n";
	OutMaterialStream << "dH0_eV=1.63;		# [eV]		enthalpy barrier for kink nucleation\n";
	OutMaterialStream << "p=0.86;				# [-]		mobility exponent\n";
	OutMaterialStream << "q=1.69;				# [-]		mobility exponent\n";
	OutMaterialStream << "Tf=0.8;				# [-]		athermal transition temperature in fraction of Tm\n";
	OutMaterialStream << "tauC_SI=2.03e9;		# [Pa]		Peierls stress\n";
	OutMaterialStream << "a0=1.50;			# [-]		non-Schmid coefficient\n";
	OutMaterialStream << "a1=1.15;			# [-]		non-Schmid coefficient\n";
	OutMaterialStream << "a2=2.32;			# [-]		non-Schmid coefficient\n";
	OutMaterialStream << "a3=4.29;			# [-]		non-Schmid coefficient\n";
	OutMaterialStream << "\n";

	OutMaterialStream << "# Vacancy diffusion\n";
	OutMaterialStream << "Omega_SI=0;	# [A^3]	Atomic volume\n";
	OutMaterialStream << "Ufv_eV=0;	# [eV]	Vacancy formation energy\n";
	OutMaterialStream << "DVv=0;	# [-]	Relative vacancy relaxation volume\n";
	OutMaterialStream << "Udv_eV=0; 	# [eV]	Vacancy migration energy\n";
	OutMaterialStream << "D0v_SI=0;	    # [m^2/s]	vacancy diffusion coefficient\n";
	*/

	
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
	
	//	Deletes the file if exists 
	if (rename(oldname, newname) != 0)
		perror("Error renaming material file\n");
	else
		std::cout << "Material file renamed successfully\n";
	
	
}

#endif
