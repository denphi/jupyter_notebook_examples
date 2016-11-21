%% OrthoCompliance.m
% Given material properties returns compliance matrix S
% Single ply of orthotropic material
% Ref. Daniel, I. M. and Ishai, O. "Engineering Mechanics of Composite
% Materials." 2nd Ed. Oxford University Press, NY (2006).
% pp 74 Eq. 4.48
%
% 10/1/07   
%
% Given a more appropreate name
% 2/12/10
%
% Andrew Ritchey
%
%%
function S = OrthoCompliance(E_1,E_2,E_3,G_12,G_13,G_23,nu_12,nu_13,nu_23)
nu_21=E_2/E_1*nu_12;
nu_31=E_3/E_1*nu_13;
nu_32=E_3/E_2*nu_23;
S=[1/E_1        -nu_21/E_2  -nu_31/E_3  0       0       0;
   -nu_12/E_1   1/E_2       -nu_32/E_3  0       0       0;
   -nu_13/E_1   -nu_23/E_2       1/E_3  0       0       0;
   0            0           0           1/G_23  0       0;
   0            0           0           0       1/G_13  0;
   0            0           0           0       0       1/G_12];