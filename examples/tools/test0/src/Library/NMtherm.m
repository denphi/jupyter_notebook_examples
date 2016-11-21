function varargout = NMtherm(Qbar,alpha,t,dT)
% Determine effective properties of a laminated plate composited of
% unidirectional composite layers.
%    Copyright (C) 2010 Andrew Ritchey
%
%    This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
%
%%  NMtherm.m
%   Calculates N_thermal and M_thermal given the Qbar, alpha, and thickness vectors for
%   the layup.  AAE 555 notes pg 85 eq. 5.63a
%   Inputs:
%       Qbar - [3x3xk] Qbar for each lamina
%       alpha - [3xk] alpha vector for each lamina
%       t - [k] thickness for each lamina
%       dT - [1x1] Temp change
%   Outputs:
%       Nt - Tensile forces due to thermal loading
%       Mt - Moments due to thermal loading

Nt=zeros(3,1);
Mt=Nt;
tt=sum(t);
z=-tt/2;

for i = 1:length(t)
    z(i+1)=z(i)+t(i);
    Nt=Nt+dT*Qbar(:,:,i)*alpha(:,i)*t(i);
    Mt=Mt+dT/2*Qbar(:,:,i)*alpha(:,i)*(z(i+1)^2-z(i)^2);
end

if nargout==1
    varargout{1}=[Nt;Mt];
elseif nargout==2
    varargout{1}=Nt;
    varargout{2}=Mt;
else
    disp('Number of output arguments should be 1 or 2.  Type "help NMtherm" for more information');
end
