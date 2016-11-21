%%  clpt.m
%   Calculates the displacement, strains and stresses in a composite 
%   laminate composed of linearly elastic lamina.
%
%   Copyright (C) 2010 Andrew Ritchey
%
%   This program is free software; you can redistribute it and/or modify 
%   it under the terms of the GNU General Public License as published by 
%   the Free Software Foundation; either version 2 of the License, or (at 
%   your option) any later version.
%
%   This program is distributed in the hope that it will be useful, 
%   but WITHOUT ANY WARRANTY; without even the implied warranty of 
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
%   General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the 
%   Free Software Foundation, Inc., 
%   59 Temple Place, Suite 330, 
%   Boston, MA 02111-1307 USA
%=======
function clpt(lib)

%path('./Library/',path)
%%  Gather Inputs
% open our xml input file.
%lib = rpLib(infile);

%%  Gather Inputs
% retrieve user specified data out of the input file
% convert values to correct units.
E1 = rpLibGetString(lib,'input.group.(Material).number(E1).current');
[E1,err] = rpUnitsConvertDbl(E1,'Pa');
E2 = rpLibGetString(lib,'input.group.(Material).number(E2).current');
[E2,err] = rpUnitsConvertDbl(E2,'Pa');

nu12 = rpLibGetString(lib,'input.group.(Material).number(nu12).current');
nu12 = str2double(nu12);

G12 = rpLibGetString(lib,'input.group.(Material).number(G12).current');
[G12,err] = rpUnitsConvertDbl(G12,'Pa');

alpha1= rpLibGetString(lib,'input.group.(Material).number(alpha1).current');
alpha1 = str2double(alpha1);

alpha2 = rpLibGetString(lib,'input.group.(Material).number(alpha2).current');
alpha2 = str2double(alpha2);

h0 = rpLibGetString(lib,'input.group.(Layup).number(h0).current');
[h0,err] = rpUnitsConvertDbl(h0,'m');

theta = rpLibGetString(lib,'input.group.(Layup).string(theta).current');
theta = str2num(theta);

Nx = rpLibGetString(lib,'input.group.(Loading).number(Nx).current');
[Nx,err] = rpUnitsConvertDbl(Nx,'N/m');

Ny = rpLibGetString(lib,'input.group.(Loading).number(Ny).current');
[Ny,err] = rpUnitsConvertDbl(Ny,'N/m');

Nxy = rpLibGetString(lib,'input.group.(Loading).number(Nxy).current');
[Nxy,err] = rpUnitsConvertDbl(Nxy,'N/m');

Mx = rpLibGetString(lib,'input.group.(Loading).number(Mx).current');
[Mx,err] = rpUnitsConvertDbl(Mx,'N');

My = rpLibGetString(lib,'input.group.(Loading).number(My).current');
[My,err] = rpUnitsConvertDbl(My,'N');

Mxy = rpLibGetString(lib,'input.group.(Loading).number(Mxy).current');
[Mxy,err] = rpUnitsConvertDbl(Mxy,'N');

NMm = [Nx;Ny;Nxy;Mx;My;Mxy];

Ti = rpLibGetString(lib,'input.group.(Loading).number(Ti).current');
[Ti,err] = rpUnitsConvertDbl(Ti,'C');

Tf = rpLibGetString(lib,'input.group.(Loading).number(Tf).current');
[Tf,err] = rpUnitsConvertDbl(Tf,'C');

dT = Tf-Ti;

%xRes = rpLibGetString(lib,'input.group.(Loading).integer(xRes).current');
%xRes = str2double(xRes);
xRes = 1;

%yRes = rpLibGetString(lib,'input.group.(Loading).integer(yRes).current');
%yRes = str2double(yRes);
yRes = 1;

[zSig,sigmaXY,sigma12,z,epsilonXY,epsilon12,surfI,surfF,laminateData] = solveCLPT(E1,E2,nu12,G12,alpha1,alpha2,h0,theta,xRes,yRes,NMm,dT);

%%  Check Inputs
outData = [E1,E2,nu12,G12,alpha1,alpha2,h0,Nx,Ny,Nxy,Mx,My,Mxy,dT];
putStr = sprintf('E1 \t %12g \n E2 \t %12g \n nu12 \t %12g \n G12 \t %12g \n alpha1 \t %12g \n alpha2 \t %12g \n h0 \t %12g \n Nx \t %12g \n Ny \t %12g \n Nxy \t %12g \n Mx \t %12g \n My \t %12g \n Mxy \t %12g \n dT \t %12g', outData);
%putStr = sprintf('%12g\n', outData);
putStr2 = sprintf('Theta %12s', num2str(theta));
rpLibPutString(lib,'output.string(Ins).current',[putStr,putStr2],0);
rpLibPutString(lib,'output.string(Ins).about.label','Inputs',0);

%%  Plot Results
%   Figure 1, Stress in laminate coordinate system through thickness
outData = [sigmaXY(1,:);zSig'/h0];
putStr = sprintf('%12g  %12g\n', outData);
rpLibPutString(lib,'output.curve(sigX).component.xy',putStr,0);
rpLibPutString(lib,'output.curve(sigX).about.group','Stress Through Thickness, Laminate Coordinate System',0);
rpLibPutString(lib,'output.curve(sigX).about.label','sigma_{x}',0);

outData = [sigmaXY(2,:);zSig'/h0];
putStr = sprintf('%12g  %12g\n', outData);
rpLibPutString(lib,'output.curve(sigY).component.xy',putStr,0);
rpLibPutString(lib,'output.curve(sigY).about.group','Stress Through Thickness, Laminate Coordinate System',0);
rpLibPutString(lib,'output.curve(sigY).about.label','sigma_{y}',0);

outData = [sigmaXY(3,:);zSig'/h0];
putStr = sprintf('%12g  %12g\n', outData);
rpLibPutString(lib,'output.curve(tauXY).component.xy',putStr,0);
rpLibPutString(lib,'output.curve(tauXY).about.group','Stress Through Thickness, Laminate Coordinate System',0);
rpLibPutString(lib,'output.curve(tauXY).about.label','tau_{xy}',0);

rpLibPutString(lib,'output.curve(tauXY).xaxis.label','Stess',0);
rpLibPutString(lib,'output.curve(tauXY).xaxis.units','Pa',0);
rpLibPutString(lib,'output.curve(tauXY).yaxis.label','Through Thickness Location z/h0',0);
rpLibPutString(lib,'output.curve(tauXY).yaxis.units','-',0);

%   Figure 2, Stress in material coordinate system through thickness
outData = [sigma12(1,:);zSig'/h0];
putStr = sprintf('%12g  %12g\n', outData);
rpLibPutString(lib,'output.curve(sig1).component.xy',putStr,0);
rpLibPutString(lib,'output.curve(sig1).about.group','Stress Through Thickness, Material Coordinate System',0);
rpLibPutString(lib,'output.curve(sig1).about.label','sigma_{1}',0);

outData = [sigma12(2,:);zSig'/h0];
putStr = sprintf('%12g  %12g\n', outData);
rpLibPutString(lib,'output.curve(sig2).component.xy',putStr,0);
rpLibPutString(lib,'output.curve(sig2).about.group','Stress Through Thickness, Material Coordinate System',0);
rpLibPutString(lib,'output.curve(sig2).about.label','sigma_{2}',0);

outData = [sigma12(3,:);zSig'/h0];
putStr = sprintf('%12g  %12g\n', outData);
rpLibPutString(lib,'output.curve(tau12).component.xy',putStr,0);
rpLibPutString(lib,'output.curve(tau12).about.group','Stress Through Thickness, Material Coordinate System',0);
rpLibPutString(lib,'output.curve(tau12).about.label','tau_{12}',0);

rpLibPutString(lib,'output.curve(tau12).xaxis.label','Stess',0);
rpLibPutString(lib,'output.curve(tau12).xaxis.units','Pa',0);
rpLibPutString(lib,'output.curve(tau12).yaxis.label','Through Thickness Location z/h0',0);
rpLibPutString(lib,'output.curve(tau12).yaxis.units','-',0);

%   Figure 3, Strain in laminate coordinate system through thickness
outData = [epsilonXY(1,:);z'/h0];
putStr = sprintf('%12g  %12g\n', outData);
rpLibPutString(lib,'output.curve(epsX).component.xy',putStr,0);
rpLibPutString(lib,'output.curve(epsX).about.group','Strain Through Thickness, Laminate Coordinate System',0);
rpLibPutString(lib,'output.curve(epsX).about.label','epsilon_{x}',0);

outData = [epsilonXY(2,:);z'/h0];
putStr = sprintf('%12g  %12g\n', outData);
rpLibPutString(lib,'output.curve(epsY).component.xy',putStr,0);
rpLibPutString(lib,'output.curve(epsY).about.group','Strain Through Thickness, Laminate Coordinate System',0);
rpLibPutString(lib,'output.curve(epsY).about.label','epsilon_{y}',0);

outData = [epsilonXY(3,:);z'/h0];
putStr = sprintf('%12g  %12g\n', outData);
rpLibPutString(lib,'output.curve(gammaXY).component.xy',putStr,0);
rpLibPutString(lib,'output.curve(gammaXY).about.group','Strain Through Thickness, Laminate Coordinate System',0);
rpLibPutString(lib,'output.curve(gammaXY).about.label','gamma_{xy}',0);

rpLibPutString(lib,'output.curve(gammaXY).xaxis.label','Strain',0);
rpLibPutString(lib,'output.curve(gammaXY).xaxis.units','-',0);
rpLibPutString(lib,'output.curve(gammaXY).yaxis.label','Through Thickness Location z/h0',0);
rpLibPutString(lib,'output.curve(gammaXY).yaxis.units','-',0);

%   Figure 4, Strain in material coordinate system through thickness
outData = [epsilon12(1,:);zSig'/h0];
putStr = sprintf('%12g  %12g\n', outData);
rpLibPutString(lib,'output.curve(eps1).component.xy',putStr,0);
rpLibPutString(lib,'output.curve(eps1).about.group','Strain Through Thickness, Material Coordinate System',0);
rpLibPutString(lib,'output.curve(eps1).about.label','epsilon_{1}',0);

outData = [epsilon12(2,:);zSig'/h0];
putStr = sprintf('%12g  %12g\n', outData);
rpLibPutString(lib,'output.curve(eps2).component.xy',putStr,0);
rpLibPutString(lib,'output.curve(eps2).about.group','Strain Through Thickness, Material Coordinate System',0);
rpLibPutString(lib,'output.curve(eps2).about.label','epsilon_{2}',0);

outData = [epsilon12(3,:);zSig'/h0];
putStr = sprintf('%12g  %12g\n', outData);
rpLibPutString(lib,'output.curve(gamma12).component.xy',putStr,0);
rpLibPutString(lib,'output.curve(gamma12).about.group','Strain Through Thickness, Material Coordinate System',0);
rpLibPutString(lib,'output.curve(gamma12).about.label','gamma_{12}',0);

rpLibPutString(lib,'output.curve(gamma12).xaxis.label','Strain',0);
rpLibPutString(lib,'output.curve(gamma12).xaxis.units','-',0);
rpLibPutString(lib,'output.curve(gamma12).yaxis.label','Through Thickness Location z/h0',0);
rpLibPutString(lib,'output.curve(gamma12).yaxis.units','-',0);
rpLibPutString(lib,'output.string(Laminate).about.label','Effective Laminate Properties',0);
%
putStr2 = sprintf('%12g \n',laminateData(1));
rpLibPutString(lib,'output.number(Ex).current',putStr2,0);
rpLibPutString(lib,'output.number(Ex).about.label','Longitudinal Modulus',0);
rpLibPutString(lib,'output.number(Ex).units','Pa',0);
%
putStr2 = sprintf('%12g \n',laminateData(2));
rpLibPutString(lib,'output.number(Ey).current',putStr2,0);
rpLibPutString(lib,'output.number(Ey).about.label','Transverse Modulus',0);
rpLibPutString(lib,'output.number(Ey).units','Pa',0);
%
putStr2 = sprintf('%12g \n',laminateData(4));
rpLibPutString(lib,'output.number(Gxy).current',putStr2,0);
rpLibPutString(lib,'output.number(Gxy).about.label','Shear Modulus',0);
rpLibPutString(lib,'output.number(Gxy).units','Pa',0);
%
putStr2 = sprintf('%12g \n',laminateData(3));
rpLibPutString(lib,'output.number(nu).current',putStr2,0);
rpLibPutString(lib,'output.number(nu).about.label','Poisson ratio',0);
%
putStr2 = sprintf('%12g \n',laminateData(5));
rpLibPutString(lib,'output.number(alphax).current',putStr2,0);
rpLibPutString(lib,'output.number(alphax).about.label','Thermal Expansion (x)',0);
%
putStr2 = sprintf('%12g \n',laminateData(6));
rpLibPutString(lib,'output.number(alphay).current',putStr2,0);
rpLibPutString(lib,'output.number(alphay).about.label','Thermal Expansion (y)',0);
%
putStr2 = sprintf('%12g \n',laminateData(7));
rpLibPutString(lib,'output.number(alphaxy).current',putStr2,0);
rpLibPutString(lib,'output.number(alphaxy).about.label','Thermal Expansion (xy)',0);

%   Effective Material Properties
outData = laminateData;
putStr2 = sprintf('Ex = %12g \n Ey = %12g \n nuxy = %12g \n Gxy = %12g \n alphax = %12g \n alphay = %12g \n alphaxy = %12g', outData);
rpLibPutString(lib,'output.string(Laminate).current',putStr2,0);

%%  Clean up
% signal the end of processing
rpLibResult(lib);
quit;
