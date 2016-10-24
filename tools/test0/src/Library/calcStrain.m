function eps = calcStrain(eps0,z)
% local eps: Strain at laminate interfaces given by z
eps=zeros(3,length(z));
for i=1:length(z)
    eps(:,i) = eps0(1:3)+z(i)*eps0(4:6);
end