function Q = ReducedStiffness(varargin)
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
%%  ReducedStiffness.m   
%   This function returns the reduced stiffness
%   matrix for fiber-reinforced materials. 
%   There are four arguments representing four  
%   material constants. The size of the reduced 
%   stiffness matrix is 3 x 3.
%   9/19/08
%   Andrew Ritchey, adapted from Matlab for composites
%   Input - Variable
%       If 1 input then must be a vector of the form 
%       [E, NU] Isotropic
%       [E1,E2,NU12,G12] (Transversly Isotropic)
%       If 2 inputs
%       E, NU (Isotropic)
%       If 4 inputs
%       E1, E2, NU12, G12 (Transversly Isotropic)
%   Outputs -
%       Q - [3x3] reduced stiffness matrix in the material coordinate
%       system.  Of the form:
%       [Q11 Q12 Q16
%        Q21 Q22 Q26
%        Q61 Q62 Q66];
%   Calls
%       Q = ReducedStiffness(isoProps);
%       Q = ReducedStiffness(laminaProps);
%       Q = ReducedStiffness(E,nu);
%       Q = ReducedStiffness(E1,E2,NU12,G12);


%%  Vector Input
if nargin==1
    % Vector of length 2, Isotropic
    if length(varargin{1})==2
        E  = varargin{1}(1);
        NU = varargin{1}(2);
        Q  = isotropic(E,NU);
    
    % Vector of length 4, Transvesly Isotropic
    elseif length(varargin{1})==4
        E1   = varargin{1}(1);
        E2   = varargin{1}(2);
        NU12 = varargin{1}(3);
        G12  = varargin{1}(4);
        Q    = transverslyIsotropic(E1,E2,NU12,G12);
    
    else
        disp('Incorrect format for input arguments, Try help ReducedStiffness for hints.')
    end

%   2 Inputs, Isotropic
elseif nargin==2
    E  = varargin{1}(1);
    NU = varargin{1}(2);
    Q  = isotropic(E,NU);

%   4 Inputs, Transversly Isotropic
elseif nargin==4
    E1   = varargin{1};
    E2   = varargin{2};
    NU12 = varargin{3};
    G12  = varargin{4};
    Q    = transverslyIsotropic(E1,E2,NU12,G12);

else
    disp('Incorrect format for input arguments, Try help ReducedStiffness for hints.')
end


function Qi = isotropic(E, NU) 
%%  Calculate the reduced stiffness matrix for an isotropic material
Qi = [E/(1-NU*NU)       NU*E/(1-NU*NU)  0 
      NU*E/(1-NU*NU)    E/(1-NU*NU)     0 
      0                 0               E/2/(1+NU)];
  
function Qti = transverslyIsotropic(E1,E2,NU12,G12)
%%  Calculate the reduced stiffness matrix for a transversly isotropic material
NU21 = NU12*E2/E1;
Qti = [E1/(1-NU12*NU21)       NU12*E2/(1-NU12*NU21)   0 
         NU12*E2/(1-NU12*NU21)  E2/(1-NU12*NU21)        0       
         0                      0                       G12];      
