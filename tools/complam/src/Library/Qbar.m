function y = Qbar(Q,theta)
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
%%  Qbar.m   
%   This function returns the transformed reduced stiffness
%   Inputs:
%       Q - [3x3] Reduced stiffness matrix in material coordinate system
%       theta - [nx1] Off axis angles in degrees for each lamina
%   Ouputs:
%       y - [3x3xn] Reduced stiffness matrix in global system for each
%       lamina
y=zeros(3,3,length(theta)); % Memory allocation
m = cos(theta*pi/180);   % cosine
n = sin(theta*pi/180);   % sine
for i=1:length(theta)       % Loop through each lamina
    Tinv = [m(i)*m(i)   n(i)*n(i)   -2*m(i)*n(i) 
            n(i)*n(i)   m(i)*m(i)   2*m(i)*n(i) 
            m(i)*n(i)   -m(i)*n(i)  m(i)*m(i)-n(i)*n(i)]; % Inverse transformation matrix
    y(:,:,i) = Tinv*Q*Tinv.';  % Qbar matrix
end
%%  Notes
%   Correct Formula, A. Ritchey 8/26/08 
%   9/16/08 Changed to a for loop
