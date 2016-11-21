function Tsig = TransSigFull(t)
%%  TransSigFull.m
%   Full 3D transformation matrix for a rotation about the 3 axis
%   5/14/08
%   Andrew Ritchey
%   Transformation matrix defined by {sigma_12}=[TSig]*{sigma_xy}
%   Inputs:
%       - Rotation angle [deg], t, 1x1
%   Outputs:
%       - Transformation matrix, [TSig], 6x6
%   Reference:
%       Sun, C. T. "Mechanics of Composite Materials and Laminates Lecture
%       Notes for AAE 555." Spring 2008.
c=cos(t*pi/180);
s=sin(t*pi/180);
Tsig=[c^2       s^2     0       0       0        2*c*s;
      s^2       c^2     0       0       0        -2*c*s;
      0         0       1       0       0        0;
      0         0       0       c       -s       0;
      0         0       0       s       c        0;
      -c*s      c*s     0       0       0        c^2-s^2;];