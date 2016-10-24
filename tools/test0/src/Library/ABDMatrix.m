function [varargout] = ABDMatrix(Qbar,z)
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
%%  ABDMatrix.m
%   Calculates the ABD matrix given the Qbar and z for each lamina
%   8/26/08
%   Andrew Ritchey
%   Inputs:
%       Qbar - [3,3,n] Qbar matrix for each lamina
%       z - [n+1] z locations (bottom at top) for each lamina


%%  Check that each var is properly sized
if size(Qbar,3)+1~=length(z)
    disp('Input size error.  See "help ABDMatrix" for more information')
end
%%  Loop over each lamina
A=zeros(3);
B=A;
D=A;
for i=1:size(Qbar,3)
    A=A+Qbar(:,:,i)*(z(i+1)-z(i));
    B=B+.5*Qbar(:,:,i)*(z(i+1)^2-z(i)^2);
    D=D+1/3*Qbar(:,:,i)*(z(i+1)^3-z(i)^3);
end
ABD=[A,B;B,D];
%%  Output
if nargout==1
    varargout{1}=ABD;
elseif nargout==3
    varargout{1}=A;
    varargout{2}=B;
    varargout{3}=D;
elseif narargout==4;
    varargout{1}=A;
    varargout{2}=B;
    varargout{3}=D;
    varargout{4}=ABD;
else
    disp('Number of output arguments should be 1, 3, or 4');
end
