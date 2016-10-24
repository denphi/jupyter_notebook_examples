function TEps = TransEps(t)
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
%% TransEps.m
%  5/14/08
%  Andrew Ritchey
%  Transformation matrix defined by {epsilon_12}=[TEps]*{epsilon_xy}
%  Inputs:
%    - Rotation angle [deg], t, 1x1
%  Outputs:
%    - Transformation matrix, [TEps], 3x3
%  Reference:
%   Sun, C. T. "Mechanics of Composite Materials and Laminates Lecture
%   Notes for AAE 555." Spring 2008.
c=cos(t*pi/180);
s=sin(t*pi/180);
TEps=[c.^2  s.^2 s.*c;
      s.^2  c.^2 -s.*c;
      -2*s.*c 2*s.*c c.^2-s.^2];
