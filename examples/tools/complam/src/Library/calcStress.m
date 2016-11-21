function [zLoc,sig] = calcStress(QBar,eps,z,alpha,dT)
%   Initialize interface locations and s
if length(z)==2 % Single Ply
    sig = zeros(3,length(z));
    zLoc = zeros(length(z),1);
else % Laminate
    sig = zeros(3,length(z)+1);
    zLoc = zeros(length(z)+1,1);
end

% Loop over plys
for i=1:size(QBar,3)
    % Strain at the bottom of the ply
    eps1 = eps(:,i);
    % Strain at the top of the ply
    eps2 = eps(:,i+1);

    % Stress in the current ply at the bottom interface
    sig(:,2*i-1) = QBar(:,:,i)*(eps1-alpha(:,i)*dT);
    zLoc(2*i-1) = z(i);

    % Stress in the current ply at the top interface
    sig(:,2*i) = QBar(:,:,i)*(eps2-alpha(:,i)*dT);
    zLoc(2*i) = z(i+1);
end