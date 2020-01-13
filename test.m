clear

dArray = [10^12];
d = 10^12;
dt = 10^-3;

a1 = (10^12)^0.1; 
a2 = (10^12)^-1.5;

mat = [0];

for n = 1:1000
    
    drate = a1*d^(0.5) - a2*d^2;
    
    d = d + drate*dt;
    
    dArray = [dArray, d];
    
    
    mat = [mat,n];
end

plot(dArray);