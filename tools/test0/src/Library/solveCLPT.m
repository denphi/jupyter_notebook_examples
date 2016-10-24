function  [zSig,sigmaXY,sigma12,z,epsilonXY,epsilon12,surfI,surfF,laminateData] = solveCLPT(E1,E2,nu12,G12,alpha1,alpha2,h0,theta,xRes,yRes,NMm,dT)
    % Q: Lamina Q matrix in material coordinate system
    % Qbar: Lamina Q matrix in laminate coordinate system
    % n: Number of plys
    n = length(theta);
    Q = zeros(3,3,n);
    QBar = Q;

    % alpha: Lamina CTE vector in lamina coordinate system
    % alphaBar: Lamina CTE vector in the laminate coordinate system
    alpha = [alpha1,alpha2,0]';
    alphaBar = zeros(3,n);

    % t: Ply by ply thickness
    % h: Total laminate thickness
    % z: Interface z-coordinates
    h0 = h0*ones(length(theta),1);
    h = sum(h0);
    z = locations(h0);
    % Loop over each ply
    for i=1:n
        % Calculate the Q matrix in the material cooridant system
        Q(:,:,i) = ReducedStiffness([E1,E2,nu12,G12]);
        % Transform the Q matrix from the material to the laminate
        % coordinate system
        QBar(:,:,i) = Qbar(Q(:,:,i),theta(i));
        % Transform the alpha vector from the material to the laminate
        % coordinate system
        alphaBar(:,i) = TransEps(-theta(i))*alpha;
    end

    % ABD: The equivalent plate stiffness matrix
    ABD = ABDMatrix(QBar,z);

    % NMt1: Induced thermal loads per unit temperature change
    NMt1 = NMtherm(QBar,alphaBar,h0,1);

    % NMt: Induced thermal loads for applied temperature change
    NMt = NMt1*dT;

    % NM: Total loads applied to laminate
    NM = NMm+NMt;

    % mechData: Mechanical laminate properties [Ex, Ey, nuxy, Gxy]
    mechData = LaminateProperties(ABD(1:3,1:3),h);

    % thermData: Thermal laminate properties [alphax,alphay,alphaxy]
    thermData = ABD\NMt1;

    % laminateData: Equivalent laminate properties [Ex, Ey, nuxy, Gxy]
    laminateData = [mechData';thermData(1:3)];

    % midStrain: Strains and curvatures at the midplane of the laminate
    % [epsilonx0, epsilony0, gammaxy0, kappax, kappay, kappaxy]
    midStrain = ABD\NM;

    % eps: Through thickness strains in the laminate
    epsilonXY = calcStrain(midStrain,z);

    
    % zPlot: Interface locations with repeats to plot the stress in each
    % ply
    % sig: Through thickness stresses in the laminate
    [zSig,sigmaXY] = calcStress(QBar,epsilonXY,z,alphaBar,dT);
    sigma12 = zeros(size(sigmaXY));
    epsilon12 = zeros(size(epsilonXY));
    
    for i=1:n
        % Transform the stress vector from the laminate to the material 
        % coordinate system
        sigma12(:,2*i-1) = TransEps(-theta(i)).'*sigmaXY(:,2*i-1);
        sigma12(:,2*i) = TransEps(-theta(i)).'*sigmaXY(:,2*i);
        
        % Determine the strain in the material coordinate system from the
        % stress in the material coordinate system.
        epsilon12(:,2*i-1) = Q(:,:,i)\sigma12(:,2*i-1);
        epsilon12(:,2*i) = Q(:,:,i)\sigma12(:,2*i);
    end

    % [xi,yi,zi]: Initial locations
    % [xf,yf,zf]: Final locations
    [surfI(:,:,1),surfI(:,:,2),surfI(:,:,3),surfF(:,:,1),surfF(:,:,2),surfF(:,:,3)] = calcDisp(midStrain,xRes,yRes);
 