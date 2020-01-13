%% Script to simulate a basic stress-strain curve

%%Simulation Constants
b = 0.5*10^-10; %Burgers constant meters
E = 30*10^9; %elastic modulus in Pa
B = 200; %arbitrary constant for dislocation velocity
strainRate = 10^-3; 
disDensity = 10^12; %starting dislocation density, in m^2
SchmidtFactor = 0.5; %Schmidt Factor/Taylor FActor
dt = 10^-3; %step magnitude for simulation

%Dislocation Density Rate Constants
a1 = (10^12)^0.5; 
a2 = (10^12)^-1.9;

%Strengthening Mechanism Constants - all given in Pa
CRSS = 1.5*10^8; %Critical Resolved Shear Stress
tauWH =0;
tauPRE =0;
tauSS = 0;
tauGB = 0;

%Vectors to hold stress and strain values
stress = [0];
strain = [0];
plasticStrain = [0];


t = 0; %time
disDensityRate = 0; %dislocation density rate of change
v = 0; %dislocation velocity

isPlastic = 0; %Binary integer to indicate when we have moved to the plastic regime

yieldStress=0;
yieldStrain=0;

for n = 1:100000
    
    t = t + dt;
    strain = [strain, strainRate*t]; %update strain vector with total strain
    
    if isPlastic ==0

        stressVal = E*strainRate*t;

        RSS = stressVal*SchmidtFactor; %Resolved shear stress
        
        totalStrengthening = CRSS + tauWH + tauPRE + tauSS + tauGB;
    
        if RSS > (totalStrengthening)
            isPlastic=1; 
            yieldStress = stressVal;
            yieldStrain = strainRate*t;
        end
        
         if RSS < (totalStrengthening)
             plasticStrain = [plasticStrain, 0];
         end    
    end
    
    
    if isPlastic ==1     
       stressVal = E*((strainRate*t + disDensity*(b^2)*totalStrengthening/B)/(1+E*disDensity*(b^2)*SchmidtFactor/B));
       
       %plasticStrainRate = (disDensity*(b^2)/B)*(RSS - totalStrengthening);
       
       disDensityRate = a1*disDensity^(0.5) - a2*disDensity^(2); %Re-calculate the dislocation density rate of change
       disDensity = disDensity + disDensityRate*dt; %Recalculate the dislocation density
    end
    
    stress = [stress,stressVal]; %update stress vector
    
    
end

scatter(strain, stress)




