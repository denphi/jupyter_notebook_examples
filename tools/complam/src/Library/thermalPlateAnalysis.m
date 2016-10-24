% %%  thermalPlateAnalysis.m
% %   10/07/08
% %   Analizes composite plates under thermal loading to predict curvature
% %   strain due to curing.
% clear 
% clc
% close all
% %%  Import FE Data
% files={'Plate2.rpt','Plate1.rpt','Plate3.rpt','Plate4.rpt'};
% layup={'[0_4/90_4]','[0/45/90/-45]','[45_2/-45_2]','[45/-45/0/90]_s'};
% for i=1:length(files)
%     importfile(files{i});
%     
%     z{:,i}=data(:,1);
%     eps_x{:,i}=data(:,2);
%     eps_y{:,i}=data(:,4);
%     eps_xy{:,i}=data(:,3);
%     
%     % Calculate Centerline Strain  
%     eps_x0{:,i}=ones(size(eps_x{:,i}))*mean(eps_x{:,i});
%     eps_y0{:,i}=ones(size(eps_x{:,i}))*mean(eps_y{:,i});
%     eps_xy0{:,i}=ones(size(eps_x{:,i}))*mean(eps_xy{:,i});
% 
%     % Calculate Curvature Induced Strain
%     eps_xk{:,i}=eps_x{:,i}-eps_x0{:,i};
%     eps_yk{:,i}=eps_y{:,i}-eps_y0{:,i};
%     eps_xyk{:,i}=eps_xy{:,i}-eps_xy0{:,i};
%     
%     % Remove insignificant figures
%     eps_x{:,i}=round(eps_x{:,i}*10^7)/10^7;
%     eps_y{:,i}=round(eps_y{:,i}*10^7)/10^7;
%     eps_xy{:,i}=round(eps_xy{:,i}*10^7)/10^7;
% 
%     eps_x0{:,i}=round(eps_x0{:,i}*10^7)/10^7;
%     eps_y0{:,i}=round(eps_y0{:,i}*10^7)/10^7;
%     eps_xy0{:,i}=round(eps_xy0{:,i}*10^7)/10^7;
% 
%     eps_xk{:,i}=round(eps_xk{:,i}*10^7)/10^7;
%     eps_yk{:,i}=round(eps_yk{:,i}*10^7)/10^7;
%     eps_xyk{:,i}=round(eps_xyk{:,i}*10^7)/10^7;
% end
% 
% clear data textdata
% 
% %%  Plot Strain Distributions
% for i=1:length(layup)
%     figure(i)
%     subplot(1,3,1)
%     plot(eps_x{:,i},z{:,i}-mean(z{:,i}),eps_y{:,i},z{:,i}-mean(z{:,i}),eps_xy{:,i},z{:,i}-mean(z{:,i}))
%     xlabel('Strain'),ylabel('Z Distance (m)')
%     title(['Strain State For ',layup{i}]);
%     subplot(1,3,2)
%     plot(eps_x0{:,i},z{:,i}-mean(z{:,i}),eps_y0{:,i},z{:,i}-mean(z{:,i}),eps_xy0{:,i},z{:,i}-mean(z{:,i}))
%     xlabel('Strain'),ylabel('Z Distance (m)')
%     title('Midplane Strain \epsilon_0');
%     subplot(1,3,3)
%     plot(eps_xk{:,i},z{:,i}-mean(z{:,i}),eps_yk{:,i},z{:,i}-mean(z{:,i}),eps_xyk{:,i},z{:,i}-mean(z{:,i}))
%     xlabel('Strain'),ylabel('Z Distance (m)')
%     legend('\epsilon_x','\epsilon_y','\epsilon_{xy}','Location','NorthEast');
%     title('Curvature Strain \epsilon_{\kappa}');
% end
% 
% %%  Print FE Results 
% t1=sprintf('%90s','          Summary of Finite Element Predicted Curvature Strains            '); 
% t2=sprintf('%20s %25s %25s %25s','Layup','\epsilon_{x\kappa} (Max)','\epsilon_{y\kappa} (Max)','\epsilon_{xy\kappa} (Max)');
% t3=sprintf('%20s %25.3e %25.3e %25.3e \n',layup{1},max(eps_xk{:,1}),max(eps_yk{:,1}),max(eps_xyk{:,1}),...
%                                           layup{2},max(eps_xk{:,2}),max(eps_yk{:,2}),max(eps_xyk{:,2}),...
%                                           layup{3},max(eps_xk{:,3}),max(eps_yk{:,3}),max(eps_xyk{:,3}),...
%                                           layup{4},max(eps_xk{:,4}),max(eps_yk{:,4}),max(eps_xyk{:,4}));
% disp(t1)
% disp(t2)
% disp(t3)

%%  Plate Theory Calculations
%   Material System 
%   Ordering: [E1; E2; E3; G12; G23; G13; nu12; nu23; nu13;
%   alpha1; alpha2; alpha3]
props_c=[118e9,7.81e9,7.81e9,3.23e9,2.50e9,3.23e9,0.323,0.562,0.323,-0.55e-6,24.3e-6,24.3e-6];

%%  Layup
theta={[45,45,-45,-45]};           % [deg] Lamina angles
   
for i=1:size(theta,1)
    t=3.175e-4*ones(size(theta{i,:}));                   % [m]   Lamina thickness

    %%  Thermal Loading
    dT=-106;                                        % [deg C]Temperature change

    %%  Calculations   %%
    h=sum(t);                                       % [m] Laminate height
    z=zeros(1,length(t)+1);                         % Memory allocation  
    z(1)=-h/2;
    for j=1:length(t)
        z(j+1)=z(j)+t(j);                           % [m]   Interface locations
    end
    
    zp{:,i}=z;
    
    Q12 = ReducedStiffness(props_c([1,2,7,4]));     % Reduced Stiffness in material coords
    QXY = Qbar(Q12,theta{i,:});                     % Reduced Stiffness for each lamina in global coords
    alphaXY = TransAlpha(...
        [props_c(10);props_c(11);0],theta{i,:});         % Coefficient of thermal expansion for each lamina in global coords.
    ABD = ABDMatrix(QXY,z);                         % Matrix which relates strain and force for a laminate
    [Ft]=NMtherm(QXY,alphaXY,t,dT);                 % Thermal in plane and bending loads under 1 unit of temperature change

    %%  Visualization on a flat plate
    epsT=inv(ABD)*Ft;                               % Calculate thermal strains for the plate
    epsT=round(epsT*10^6)/10^6;                     % Remove numerical error
    
    eps_x0A{:,i}=ones(size(z)).'*epsT(1);
    eps_y0A{:,i}=ones(size(z)).'*epsT(2);
    eps_xy0A{:,i}=ones(size(z)).'*epsT(3);
    
    eps_xkA{:,i}=z.'*epsT(4);
    eps_ykA{:,i}=z.'*epsT(5);
    eps_xykA{:,i}=z.'*epsT(6);
    
    eps_xA{:,i}=eps_x0A{:,i}+eps_xkA{:,i};
    eps_yA{:,i}=eps_y0A{:,i}+eps_ykA{:,i};
    eps_xyA{:,i}=eps_xy0A{:,i}+eps_xykA{:,i};
end
[x,y]=meshgrid(-0.2:0.01:0.2);       % Make a grid in the xy plane

scale=1;                                        % Displacement scale factor
z=zeros(size(x));                               % Initial z-height is 0
dx=epsT(1)*x+epsT(3)/2*y;                       % Change in x direction
dy=epsT(2)*y+epsT(3)/2*x;                       % Change in y direction
dz=-epsT(4)*x.^2/2-epsT(5)*y.^2/2-2*epsT(6)*x.*y;     % Change in z direction
xf=x+dx;                                        % Final x locations
yf=y+dy;                                        % Final y locations
zf=z+dz;                                        % Final z locations
figure(1)
mesh(x,y,z),hold all                            % Plot initial plate
surf(xf,yf,zf),axis equal;                      % Plot final plate
legend('Initial','Final')
disp(max(max(zf)));


[x,y]=meshgrid(0.0:0.01:0.4);       % Make a grid in the xy plane
scale=1;                                        % Displacement scale factor
z=zeros(size(x));                               % Initial z-height is 0
dx=epsT(1)*x+epsT(3)/2*y;                       % Change in x direction
dy=epsT(2)*y+epsT(3)/2*x;                       % Change in y direction
dz=-epsT(4)*x.^2/2-epsT(5)*y.^2/2-2*epsT(6)*x.*y;     % Change in z direction
xf=x+dx;                                        % Final x locations
yf=y+dy;                                        % Final y locations
zf=z+dz;                                        % Final z locations
figure(2)
mesh(x,y,z),hold all                            % Plot initial plate
surf(xf,yf,zf),axis equal;                      % Plot final plate
legend('Initial','Final')
disp(max(max(zf)));

[x,y]=meshgrid(-0.4:0.01:0.4);       % Make a grid in the xy plane
scale=1;                                        % Displacement scale factor
z=zeros(size(x));                               % Initial z-height is 0
dx=epsT(1)*x+epsT(3)/2*y;                       % Change in x direction
dy=epsT(2)*y+epsT(3)/2*x;                       % Change in y direction
dz=-epsT(4)*x.^2/2-epsT(5)*y.^2/2-2*epsT(6)*x.*y;     % Change in z direction
xf=x+dx;                                        % Final x locations
yf=y+dy;                                        % Final y locations
zf=z+dz;                                        % Final z locations
figure(3)
mesh(x,y,z),hold all                            % Plot initial plate
surf(xf,yf,zf),axis equal;                      % Plot final plate
legend('Initial','Final')
disp(max(max(zf)));

% %%  Plot Strain Distributions
% for i=1:length(layup)
%     figure(8+i)
%     subplot(1,3,1)
%     plot(eps_xA{:,i},zp{:,i},eps_yA{:,i},zp{:,i},eps_xyA{:,i},zp{:,i})
%     xlabel('Strain'),ylabel('Z Distance (m)')
%     title(['Strain State For ',layup{i}]);
%     subplot(1,3,2)
%     plot(eps_x0A{:,i},zp{:,i},eps_y0A{:,i},zp{:,i},eps_xy0A{:,i},zp{:,i})
%     xlabel('Strain'),ylabel('Z Distance (m)')
%     title('Midplane Strain \epsilon_0');
%     subplot(1,3,3)
%     plot(eps_xkA{:,i},zp{:,i},eps_ykA{:,i},zp{:,i},eps_xykA{:,i},zp{:,i})
%     xlabel('Strain'),ylabel('Z Distance (m)')
%     legend('\epsilon_x','\epsilon_y','\epsilon_{xy}','Location','NorthEast');
%     title('Curvature Strain \epsilon_{\kappa}');
% end
% 
% %%  Print Analytic Results
% t1=sprintf('%90s','          Summary of Plate Theory Predicted Curvature Strains            '); 
% t2=sprintf('%20s %25s %25s %25s','Layup','\epsilon_{x\kappa} (Max)','\epsilon_{y\kappa} (Max)','\epsilon_{xy\kappa} (Max)');
% t3=sprintf('%20s %25.3e %25.3e %25.3e \n',layup{1},max(eps_xkA{:,1}),max(eps_ykA{:,1}),max(eps_xykA{:,1}),...
%                                           layup{2},max(eps_xkA{:,2}),max(eps_ykA{:,2}),max(eps_xykA{:,2}),...
%                                           layup{3},max(eps_xkA{:,3}),max(eps_ykA{:,3}),max(eps_xykA{:,3}),...
%                                           layup{4},max(eps_xkA{:,4}),max(eps_ykA{:,4}),max(eps_xykA{:,4}));
% disp(t1)
% disp(t2)
% disp(t3)