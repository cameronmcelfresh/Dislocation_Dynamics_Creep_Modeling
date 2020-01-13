clc
close all
clear all

MODEL_DIR='../../../..';

filename='cube';
meshID=1; % creates ../N/N_1.txt and ../T/T_1.txt
L=4000;     % side length of the cube
t=150;       % tickness of the boundary layer
x0=0;     % center of cube (x coordinate)
y0=0;     % center of cube (y coordinate)
z0=0;     % center of cube (z coordinate)

%% Define vertices
% Coordinates of the outer cube
P(1,:)=[0 0 0]*L ;
P(2,:)=[1 0 0]*L;
P(3,:)=[1 1 0]*L;
P(4,:)=[0 1 0]*L;
P(5,:)=[0 0 1]*L;
P(6,:)=[1 0 1]*L;
P(7,:)=[1 1 1]*L;
P(8,:)=[0 1 1]*L;

% Coordinates of the inner cube
P( 9,:)=[t   t t];
P(10,:)=[L-t  t t];
P(11,:)=[L-t L-t t];
P(12,:)=[t  L-t t];
P(13,:)=[t   t L-t];
P(14,:)=[L-t  t L-t];
P(15,:)=[L-t L-t L-t];
P(16,:)=[t  L-t L-t];

% Recenter
P=P+repmat([-L/2+x0 -L/2+y0 -L/2+z0],size(P,1),1);

figure(1)
clf
plot3(P(:,1),P(:,2),P(:,3),'ro','Linewidth',2)
hold on
text(P(:,1),P(:,2),P(:,3),num2str([1:size(P,1)]'),'FontSize',16)
axis equal
xlabel('x')
ylabel('y')
grid on
%print(gcf,'-depsc','../dox/vertices')

%% Define facets
% Facets normal to x
Fx=[ 1  5  8  4;
     2  6  7  3;
     9 13 16 12;
    10 14 15 11];
% Facets normal to y
Fy=[ 1  2  6  5;
     9 10 14 13;
    12 11 15 16;
     4  3  7  8];
% Facets normal to z
Fz=[ 1  2  3  4;
     5  6  7  8;
     9 10 11 12;
    13 14 15 16];
 % Inclined facets
Fxy=[1 9 13 5;
    2 10 14 6;
    3 11 15 7;
    4 12 16 8];
Fxz=[1 9 12 4;
    2 10 11 3;
    6 14 15 7;
    5 13 16 8];
Fyz=[1 9 10 2;
    5 13 14 6;
    7 15 16 8;
    4 12 11 3];

F=[Fx;Fy;Fz;Fxy;Fxz;Fyz];
for f=1:size(F,1)
    vID=F(f,:);
    fill3(P(vID,1),P(vID,2),P(vID,3),'g','Facealpha',0.2)
    drawnow
    pause(0.1)
end
%print(gcf,'-depsc','../dox/zFaces')

%% Define Regions
% Each row in the following matrix correponds to a region.
% Each row contains the vertices that bound that region.
RP=[9 10 11 12 13 14 15 16;
    1 2 6 5 9 10 14 13;
    2 6 7 3 10 11 14 15;
   4 3 7 8 11 12 15 16;
   1 4 5 8 9 13 16 12;
   1 2 3 4 9 10 11 12;
   5 6 7 8 13 14 15 16];

% The rows of the matrix Xm are the barycenters of each region
for r=1:size(RP,1)
Xm(r,:)=mean(P(RP(r,:),:));
end
plot3(Xm(:,1),Xm(:,2),Xm(:,3),'bx','Linewidth',2)
text(Xm(:,1),Xm(:,2),Xm(:,3),num2str([1:size(Xm,1)]'),'FontSize',16)

%% Write .poly file
polyFile = fopen([filename '.poly'],'w');

% Part 1- the node list.
pointFormat='%i %d %d %d \n';
fprintf(polyFile,'# Part 1 - the node list.\n');
fprintf(polyFile,'%i 3 0 0 \n',size(P,1));  % number of nodes
for k=1:size(P,1)
fprintf(polyFile,pointFormat,[k P(k,:)]);
end
fprintf(polyFile,'\n');
fprintf(polyFile,'\n');

% Part 2- the facet list.
fprintf(polyFile,'# Part 2 - the facet list.\n');
fprintf(polyFile,'%i 0 \n',size(F,1)); % number of facets
for r=1:size(F,1)
fprintf(polyFile,'1 \n'); 
fprintf(polyFile,'4   %d %d %d %d \n',F(r,:));
end

% Part 3- the hole list.
fprintf(polyFile,'# Part 3 - the hole list.\n');
fprintf(polyFile,'0 \n\n');

% Part 4- the region list.
fprintf(polyFile,'# Part 4 - the region list.\n');
fprintf(polyFile,'%i \n',size(Xm,1));

meshSize=ones(size(Xm,1),1)*10000000;
meshSize(1)=meshSize(1);
for r=1:size(Xm,1)
%fprintf(polyFile,'%i %d %d %d %d \n',[r Xm(r,:) meshSize(r)]);
fprintf(polyFile,'%i %d %d %d %i %d \n',[r Xm(r,:) r meshSize(r)]);
end

fclose(polyFile);

%% Run Tetgen 
system([MODEL_DIR '/scripts/tetgenPOLY.sh ' filename]);

%% Create T and N files and clean tetgent output
system([MODEL_DIR '/scripts/tetgen2TN.sh ' filename ' ' num2str(meshID)]);
