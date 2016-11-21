function z = locations(t)
%%  locations.m
%   Given a vector of thicknesses, return the location of the interface for
%   each ply.
%
%   2/12/10
%   Andrew Ritchey

%%
z = zeros(length(t)+1,1);
tTot = sum(t);
z(1) = -tTot/2;
for i=1:length(t)
        z(i+1) = z(i)+t(i);
end