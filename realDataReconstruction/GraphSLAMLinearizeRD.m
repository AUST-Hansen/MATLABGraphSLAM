function SLAMdataOut = GraphSlamLinearizeRD(SLAMdata,OtherData)

NmapObservations = size(OtherData.Descriptors,2);
NmapFeatures = sum(unique(SLAMdata.c_i_t)>0);
stateEstimate = reshape(SLAMdata.mu(1:SLAMdata.stateSize*SLAMdata.Nstates),SLAMdata.stateSize,[]);
mapEstimate = reshape(SLAMdata.mu(SLAMdata.stateSize*SLAMdata.Nstates+1:end),SLAMdata.mapDimension,[]);
%% Initialize Information matrix and vector
% state has position, velocity, and heading info, as well as berg omega
% estimate. Map features have 3d position.
Omega = spalloc((SLAMdata.stateSize*SLAMdata.Nstates+SLAMdata.mapDimension*NmapFeatures),(SLAMdata.stateSize*SLAMdata.Nstates+SLAMdata.mapDimension*NmapFeatures),10e7);
%Omega = eps*speye(stateSize*Nstates+3*NmapFeatures);
% anchor initial position and orientation
Omega(1,1) = 1e13;
Omega(2,2) = 1e13;
Omega(3,3) = 0;
Omega(4,4) = 0;
Omega(5,5) = 1e13;
Omega(6,6) = 1e13; % TODO: reduce this
% ininitialize zeta
zeta = zeros(SLAMdata.stateSize*SLAMdata.Nstates + SLAMdata.mapDimension*NmapFeatures,1);
zeta(1:6) = [0 0 OtherData.inputs(1,1) 0 0 0]';
xhat = zeros(SLAMdata.stateSize,1);
Gmat = zeros(SLAMdata.stateSize);
%processNoise = diag([1e-10,1e-10,1e-9,1e-9,1e-9,1e-8]);
processNoise = zeros(SLAMdata.stateSize);
processNoise(1:2,1:2) = 1e-7*eye(2);
processNoise(3,3) = 1e-6;
processNoise(4,4) = 1e-6;
processNoise(5,5) = 1e-10;     % gyro noise
processNoise(end,end) = 1e-8; % bias driver

%% Motion model
fprintf('\nMotion...\n')
for ii = 1:SLAMdata.Nstates-1
    dT = (OtherData.timeSteps(ii+1) - OtherData.timeSteps(ii));
    
    %% Extract heading
    % nonlinear state update
    %[~, xHat] = ode45(@dStateDt,[0 dT],[stateEstimate(:,ii);OtherData.inputs(:,ii)]);
    %xhat = xHat(end,1:SLAMdata.stateSize)';
    dX = dStateDt(1,[stateEstimate(:,ii);OtherData.inputs(:,ii)])*dT;
    xhat = stateEstimate(:,ii) + dX(1:SLAMdata.stateSize);
    %[xhat-stateEstimate(:,ii)]
    % current state estimate
    Gmat = buildGmatrix(stateEstimate(:,ii),dT);
    % add information to Omega and zeta
    Omega((ii-1)*SLAMdata.stateSize+1:(ii+1)*SLAMdata.stateSize,(ii-1)*SLAMdata.stateSize+1:(ii+1)*SLAMdata.stateSize) = ... %self
        Omega((ii-1)*SLAMdata.stateSize+1:(ii+1)*SLAMdata.stateSize,(ii-1)*SLAMdata.stateSize+1:(ii+1)*SLAMdata.stateSize)...
        + [-Gmat'; eye(SLAMdata.stateSize)]*(processNoise\[-Gmat, eye(SLAMdata.stateSize)]);

    zeta((ii-1)*SLAMdata.stateSize+1:(ii+1)*SLAMdata.stateSize) = zeta((ii-1)*SLAMdata.stateSize+1:(ii+1)*SLAMdata.stateSize) +...
        [-Gmat'; eye(SLAMdata.stateSize)]*(processNoise\(xhat - Gmat*stateEstimate(:,ii)));
    if (mod(ii,500)==0)
        fprintf('%d of %d points...\n',ii,SLAMdata.Nstates)
    end
end

%% Measurements

fprintf('\nMeasurements...\n')
% loop over times
for t = 1:SLAMdata.Nstates
    measurements_t = SLAMdata.c_i_t(SLAMdata.c_i_t(:,t)>0,t);
    % test to see if we saw anything this time 'round
    if ~isempty(measurements_t)
        for jj = 1:size(measurements_t,1)
        % for each measurement recorded at time t
        % build zhat
        zMeas = OtherData.Descriptors(1:4,OtherData.measIndex(jj,t));
        zhat = buildZhat(stateEstimate(:,t),mapEstimate(:,measurements_t(jj)),OtherData.Z(t));
        H = buildHmatrix(stateEstimate(:,t),mapEstimate(:,measurements_t(jj)));
        mapFeatureIterator = measurements_t(jj);
        Qinv = 100*eye(4);
        CovAdd = H'*Qinv*H;
        zetaAdd = H'*Qinv*(zMeas - zhat + H*[stateEstimate(:,t); mapEstimate(:,measurements_t(jj))]);
        % update information matrix/vector
        % traj diagonal
        Omega(t*SLAMdata.stateSize-(SLAMdata.stateSize-1):t*SLAMdata.stateSize,...
              t*SLAMdata.stateSize-(SLAMdata.stateSize-1):t*SLAMdata.stateSize) = ... %self
         Omega(t*SLAMdata.stateSize-(SLAMdata.stateSize-1):t*SLAMdata.stateSize,...
               t*SLAMdata.stateSize-(SLAMdata.stateSize-1):t*SLAMdata.stateSize)...
            + CovAdd(1:SLAMdata.stateSize,1:SLAMdata.stateSize);
        % traj map

        Omega(t*SLAMdata.stateSize-(SLAMdata.stateSize-1):t*SLAMdata.stateSize,... % pos
                SLAMdata.Nstates*SLAMdata.stateSize+1+(mapFeatureIterator-1)*(SLAMdata.mapDimension):...
                SLAMdata.Nstates*SLAMdata.stateSize+mapFeatureIterator*SLAMdata.mapDimension) = ... %self
         Omega(t*SLAMdata.stateSize-(SLAMdata.stateSize-1):t*SLAMdata.stateSize,... % pos
                SLAMdata.Nstates*SLAMdata.stateSize+1+(mapFeatureIterator-1)*(SLAMdata.mapDimension):...
                SLAMdata.Nstates*SLAMdata.stateSize+mapFeatureIterator*SLAMdata.mapDimension)...
            + CovAdd(1:SLAMdata.stateSize,SLAMdata.stateSize+1:SLAMdata.stateSize+SLAMdata.mapDimension);
        % traj zeta
        zeta(t*SLAMdata.stateSize-(SLAMdata.stateSize-1):t*SLAMdata.stateSize) =...
            zeta(t*SLAMdata.stateSize-(SLAMdata.stateSize-1):t*SLAMdata.stateSize) +...
            + zetaAdd(1:SLAMdata.stateSize);

        %          % Add information about map to Omega and zeta
        % map diagonal
        Omega(SLAMdata.Nstates*SLAMdata.stateSize+1+(mapFeatureIterator-1)*(SLAMdata.mapDimension):...
              SLAMdata.Nstates*SLAMdata.stateSize+(mapFeatureIterator)*(SLAMdata.mapDimension),...
              SLAMdata.Nstates*SLAMdata.stateSize+1+(mapFeatureIterator-1)*(SLAMdata.mapDimension):...
              SLAMdata.Nstates*SLAMdata.stateSize+(mapFeatureIterator)*(SLAMdata.mapDimension)) = ... %self
           Omega(SLAMdata.Nstates*SLAMdata.stateSize+1+(mapFeatureIterator-1)*(SLAMdata.mapDimension):...
              SLAMdata.Nstates*SLAMdata.stateSize+(mapFeatureIterator)*(SLAMdata.mapDimension),...
              SLAMdata.Nstates*SLAMdata.stateSize+1+(mapFeatureIterator-1)*(SLAMdata.mapDimension):...
              SLAMdata.Nstates*SLAMdata.stateSize+(mapFeatureIterator)*(SLAMdata.mapDimension))...
              + sparse(CovAdd(SLAMdata.stateSize+1:end,SLAMdata.stateSize+1:end));

        Omega(SLAMdata.Nstates*SLAMdata.stateSize+1+(mapFeatureIterator-1)*(SLAMdata.mapDimension):...
              SLAMdata.Nstates*SLAMdata.stateSize+mapFeatureIterator*SLAMdata.mapDimension,...
              t*SLAMdata.stateSize-(SLAMdata.stateSize-1):t*SLAMdata.stateSize) = ...
            Omega(SLAMdata.Nstates*SLAMdata.stateSize+1+(mapFeatureIterator-1)*(SLAMdata.mapDimension):...
                  SLAMdata.Nstates*SLAMdata.stateSize+mapFeatureIterator*SLAMdata.mapDimension,...
                  t*SLAMdata.stateSize-(SLAMdata.stateSize-1):t*SLAMdata.stateSize) ...
            + CovAdd(SLAMdata.stateSize+1:SLAMdata.stateSize+SLAMdata.mapDimension,1:SLAMdata.stateSize);
        
        zeta(SLAMdata.Nstates*SLAMdata.stateSize+1+(mapFeatureIterator-1)*(SLAMdata.mapDimension):...
              SLAMdata.Nstates*SLAMdata.stateSize+(mapFeatureIterator)*(SLAMdata.mapDimension)) = ...
              zeta(SLAMdata.Nstates*SLAMdata.stateSize+1+(mapFeatureIterator-1)*(SLAMdata.mapDimension):...
              SLAMdata.Nstates*SLAMdata.stateSize+(mapFeatureIterator)*(SLAMdata.mapDimension)) ...
              + zetaAdd(SLAMdata.stateSize+1:end);
        %end
        end % for jj
    end % if ~isempty
    
    QDVL = 1e-6;
    QinvDVL = eye(2)*1/QDVL;
    Hdvl = [0 0 1 0 0 0; 0 0 0 1 0 0];
    zMeas = OtherData.DVL(1:2,t);
    %zHat = Hdvl*stateEstimate(:,t);
    CovAdd = Hdvl'*QinvDVL*Hdvl;
    zetaAdd = Hdvl'*QinvDVL*(zMeas);
    Omega((t-1)*SLAMdata.stateSize+1:(t)*SLAMdata.stateSize,(t-1)*SLAMdata.stateSize+1:(t)*SLAMdata.stateSize) = ...
        Omega((t-1)*SLAMdata.stateSize+1:(t)*SLAMdata.stateSize,(t-1)*SLAMdata.stateSize+1:(t)*SLAMdata.stateSize) + ...
        CovAdd;
    zeta((t-1)*SLAMdata.stateSize+1:(t)*SLAMdata.stateSize) = zeta((t-1)*SLAMdata.stateSize+1:(t)*SLAMdata.stateSize) +...
                    + zetaAdd;
    
    
    if (mod(t,500)==0)
        fprintf('%d of %d points...\n',t,SLAMdata.Nstates)
    end
end % for t

SLAMdataOut = SLAMdata;
SLAMdataOut.Omega = Omega;
SLAMdataOut.zeta = zeta;
end

%% State-specific functions


    function deriv = dStateDt(~,stateAndInputs)
        
        head = stateAndInputs(5);
        R = Euler2RotMat(0,0,head);
        b_R_v = R(1:2,1:2);
        
        deriv = zeros(size(stateAndInputs));
        deriv(1:2) = b_R_v*stateAndInputs(3:4);
        deriv(3:4) = .001*[stateAndInputs(7); 0] - stateAndInputs(3:4);
        deriv(5) = (stateAndInputs(8) - stateAndInputs(6));
        deriv(6) = -.1*stateAndInputs(6);
        % zeros for everything else
    end

    function Gmat = buildGmatrix(state,dT)
        u = state(3);
        v = state(4);
        cPsi = cos(state(5));
        sPsi= sin(state(5));
        R = [cPsi, -sPsi; sPsi, cPsi];
        Gmat(1:2,1:2) = eye(2);
        Gmat(1:2,3:4) = R*dT;
        Gmat(1:2,5) = dT*([-sPsi, -cPsi; cPsi, -sPsi]*[u; v]);
        Gmat(3:4,3:4) = eye(2);
        Gmat(5:6,5:6) = .9*eye(2);
        Gmat(5,6) = -dT;
    end

    function zHat = buildZhat(stateEstimate,mapEstimate,Z)
        
        head = stateEstimate(5);
        iRv = Euler2RotMat(0,0,head);
        x = [stateEstimate(1:2);Z];
        mapLoc = mapEstimate(1:3);
        
        zHat = [iRv'*(mapLoc-x); mapEstimate(4)];
         
        
    end

    function Hmat = buildHmatrix(state,map)
        
        head = state(5);
        iRv = Euler2RotMat(0,0,head);
        vRi = iRv';
        dmeasdhead = [-sin(head) cos(head);-cos(head) -sin(head)]*(map(1:2)-state(1:2));
        A = [ [-vRi(1:2,1:2), zeros(2,2), dmeasdhead zeros(2,1); zeros(1,6) ], iRv'];
        Hmat = blkdiag(A,1);
        
    end