clear
clc
close all

%	Material Inputs
E1 = 150e9;			% [Pa]
E2 = 9e9;			% [Pa]
nu12 = 0.33;			% [-]
G12 = 3.88e9;			% [Pa]
alpha1 = -0.9e-6;		% [strain/deg C]
alpha2 = 10e-6;			% [strain/deg C]
%theta = [45; -45; 45; -45; -45; 45; -45; 45];		% [deg]
theta = [45; -45; -45; 45];		% [deg]
t = 3.94e-3/4*ones(size(theta));	% [m]

%	Loading Inputs
dT = 0;		% [deg C]
Nm = [1;0;0];		% [N/m]
Mm = [0;0;0];		% [N]

%	Initialize variables
plyData = zeros(length(theta),8);
Q = zeros(3,3,size(plyData,1));
QBarV = Q;
TSig = Q;
alpha = zeros(3,size(plyData,1));

%	Build Laminate
for i=1:length(theta)
	plyData(i,:) = [E1,E2,nu12,G12,alpha1,alpha2,t(i),theta(i)];
end

%	Calculations
for i=1:size(plyData,1)
	Q(:,:,i) = ReducedStiffness(plyData(i,1:4));
	QBarV(:,:,i) = Qbar(Q(:,:,i),plyData(i,end));
	alpha(:,i) = TransEps(-plyData(i,end))*[plyData(i,5:6),0]';
    TSig(:,:,i) = TransEps(-plyData(i,end))';
end

h = sum(t);				% [m] Total thickness
z = locations(t);			% [m] Interface locations
ABD = ABDMatrix(QBarV,z);		% Force/Strain Matrix		
[Nt,Mt] = NMtherm(QBarV,alpha,t,dT);	% Thermal edge force/moments 
epsKap = ABD\[Nm+Nt;Mm+Mt];		% Curvature/Strain

eps0 = epsKap(1:3);			% Midline Strain
kap = epsKap(4:6);			% Curvature

epsTot = calcStrain(epsKap,z);	% Total strain through thickness
[zSig,sig] = calcStress(QBarV,epsTot,z,alpha,dT);	% Stress through 

sig12 = zeros(size(sig));
for i=1:size(sig,2)
    sig12(:,i) = TSig(:,:,ceil(i/2))*sig(:,i);
end
% [xi,yi,zi,xf,yf,zf] = calcDisp(eps0,kap);		% Displacments
 
%	Symmetric Layup Laminate Properties
lamMProps = LaminateProperties(ABD(1:3,1:3),h);
lamTProps = ABD(1:3,1:3)\(Nt/dT);

%	Output
outData = [zSig';sig(1,:);sig(2,:);sig(3,:);sig12(1,:);sig12(2,:);sig12(3,:);];
tStr = sprintf('%12s %12s %12s %12s %12s %12s %12s','z','sig_x','sig_y','tau_xy','sig_1','sig_2','tau_12'); 
putStr = sprintf('%12g %12g %12g %12g %12g %12g %12g \n',outData);
disp(tStr)
disp(putStr)