function varargout=LaminateProperties(A,h)
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
%%  LaminateProperties.m
%   4/24/08
%   Andrew Ritchey
%   Laminate properties for balanced symmetric 
%   laminates.
%   Inputs:
%       A - A matrix
%       h - Total thickness of laminate
%   Outputs:
%       Ex - Young's modulus in material x
%       Ey - Young's modulus in material y
%       Gxy - Shear stiffness in xy
%       nuxy - Poissons ratio loaded in strain y/strain x
%   Can be output as a single row vector or four individual variables.
a = inv(A);
Ex = 1/(h*a(1,1));
Ey = 1/(h*a(2,2));
Gxy = 1/(h*a(3,3));
nuxy = -a(1,2)/a(1,1);
if nargout==1,
    varargout={[Ex,Ey,nuxy,Gxy]};
elseif nargout==4,
    varargout(1)=Ex;
    varargout(2)=Ey;
    varargout(3)=nuxy;
    varargout(4)=Gxy;
end
