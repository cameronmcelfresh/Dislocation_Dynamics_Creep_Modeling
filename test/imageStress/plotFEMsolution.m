clc
clear all
close all

step=90;
U=load(['U/U_' num2str(step) '.txt']);
S=load(['S/S_' num2str(step) '.txt']);

%% Plot displacement 
comp={'x_1','x_2','x_3','u_1','u_2','u_3'};
figure(1)
hold on
axis equal

clrCol=4;
uMax=max(U(:,clrCol))
uMin=min(U(:,clrCol))
clrID=round((U(:,clrCol)-uMin)/(uMax-uMin)*size(colormap,1));

def=100;

for f=1:3:size(S,1)
fill3(U(f:f+2,1)+def*U(f:f+2,4),U(f:f+2,2)+def*U(f:f+2,5),U(f:f+2,3)+def*U(f:f+2,6),clrID(f:f+2))
plot3(U(f:f+2,1),U(f:f+2,2),U(f:f+2,3),'k')
end
grid on
xlabel('X_1','FontSize',14)
ylabel('X_2','FontSize',14)
zlabel('X_3','FontSize',14)
title(comp{clrCol},'FontSize',14)

%% Plot stress 
comp={'x_1','x_2','x_3','\sigma_{11}','\sigma_{22}','\sigma_{33}','\sigma_{12}','\sigma_{23}','\sigma_{13}'};
figure(2)
hold on
axis equal

clrCol=9;
sMax=max(S(:,clrCol))
sMin=min(S(:,clrCol))
clrID=round((S(:,clrCol)-sMin)/(sMax-sMin)*size(colormap,1));

for f=1:3:size(S,1)
fill3(S(f:f+2,1),S(f:f+2,2),S(f:f+2,3),clrID(f:f+2))
end
grid on
xlabel('X_1','FontSize',14)
ylabel('X_2','FontSize',14)
zlabel('X_3','FontSize',14)
title(comp{clrCol},'FontSize',14)

