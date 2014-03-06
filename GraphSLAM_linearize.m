function [Omega,zeta,c_i_t_new] = GraphSLAM_linearize(timeStamps,inputs,measTimestamps,sensor,rangeMeasurements,dvl,dvlReadings,correspondences,meas_ind,fullStateEstimate,haveInitialized,useRelHeading,relHeading,LCobj)

Nstates = size(timeStamps,2);
stateSize = 6;
NmapObservations = size(rangeMeasurements,2);
NmapFeatures = sum(unique(correspondences)>0); %max(max(correspondences)); %NmapObservations;
stateEstimate = reshape(fullStateEstimate(1:stateSize*Nstates),stateSize,[]);
mapEstimate = reshape(fullStateEstimate(stateSize*Nstates+1:end),3,[]);
% state has position, velocity, and heading info, as well as berg omega
% estimate. Map features have 3d position.
Omega = spalloc((stateSize*Nstates+3*NmapFeatures),(stateSize*Nstates+3*NmapFeatures),10e7);
%Omega = eps*speye(stateSize*Nstates+3*NmapFeatures);
% anchor initial position and orientation
Omega(1,1) = 1e13;
Omega(2,2) = 1e13;
Omega(3,3) = 0;
Omega(4,4) = 0;
Omega(5,5) = 1e13;
Omega(6,6) = 1; % TODO: reduce this
% ininitialize zeta
zeta = zeros(stateSize*Nstates + 3*NmapFeatures,1);
zeta(1:6) = [0 0 inputs(1,1) 0 0 0]';
xhat = zeros(stateSize,1);
G = zeros(stateSize);
%processNoise = diag([1e-10,1e-10,1e-9,1e-9,1e-9,1e-8]);
processNoise = zeros(stateSize);
processNoise(1:2,1:2) = 1e-12*eye(2);
processNoise(3,3) = 1e-4;
processNoise(4,4) = 1e-4;
processNoise(5,5) = 1e-8;     % gyro noise
processNoise(end,end) = 1e-8; % bias driver
QinvSM =5;
%% Motion model
for ii = 1:Nstates-1
    dT = (timeStamps(ii+1) - timeStamps(ii));
    heading = stateEstimate(5,ii);
    berg_R_veh = [cos(heading), -sin(heading); sin(heading), cos(heading)];

    %xhat(1:2) = eye(2)*stateEstimate(1:2,ii) + berg_R_veh*stateEstimate(3:4,ii)*dT; % update position
    %xhat(3:4) = .0*stateEstimate(3:4,ii) + 1*[inputs(1,ii); 0];
    %xhat(5) = stateEstimate(5,ii) + (inputs(2,ii) - stateEstimate(6,ii))*dT;
    %xhat(6) = stateEstimate(6,ii);
    
    [~, xHat] = ode45(@dStateDt,[0 dT],[stateEstimate(:,ii);inputs(:,ii)]);
    xhat = xHat(end,1:6)';
    % current state estimate
    u = stateEstimate(3,ii);
    v = stateEstimate(4,ii);
    bias = stateEstimate(6,ii);
    cPsi = cos(heading);
    sPsi= sin(heading);
    G(1:2,1:2) = eye(2);
    G(1:2,3:4) = berg_R_veh*dT;
    G(1:2,5) = dT*([-sPsi, -cPsi; cPsi, -sPsi]*[u; v]);
    G(3:4,3:4) = eye(2);
    G(5:6,5:6) = eye(2);
    G(5,6) = -dT;
    % add information to Omega and zeta
    Omega((ii-1)*stateSize+1:(ii+1)*stateSize,(ii-1)*stateSize+1:(ii+1)*stateSize) = ... %self
        Omega((ii-1)*stateSize+1:(ii+1)*stateSize,(ii-1)*stateSize+1:(ii+1)*stateSize)...
        + [-G'; eye(stateSize)]*(processNoise\[-G, eye(stateSize)]);

    zeta((ii-1)*stateSize+1:(ii+1)*stateSize) = zeta((ii-1)*stateSize+1:(ii+1)*stateSize) +...
        [-G'; eye(stateSize)]*(processNoise\(xhat - G*stateEstimate(:,ii)));
    
    if(useRelHeading)
        fprintf('relHeading\n')
        for jj = ii+1:Nstates
            if (relHeading.validMatch(ii,jj) == 1)
                if (abs(stateEstimate(5,jj) - stateEstimate(5,ii) - 2*pi) < .5 ) % deal with wrapping
                    % diags
                    Omega(stateSize*ii-1,stateSize*ii-1) = Omega(stateSize*ii-1:stateSize*ii-1) + QinvSM;
                    Omega(stateSize*jj-1,stateSize*jj-1) = Omega(stateSize*jj-1:stateSize*jj-1) + QinvSM;
                    % offdiags
                    Omega(stateSize*ii-1,stateSize*jj-1) = Omega(stateSize*ii-1:stateSize*ii-1) - QinvSM;
                    Omega(stateSize*jj-1,stateSize*ii-1) = Omega(stateSize*jj-1:stateSize*jj-1) - QinvSM;
                    % zeta
                    zeta(stateSize*jj-1) = zeta(stateSize*jj-1) + QinvSM*(relHeading.deltaHeading(ii,jj)+2*pi);
                    zeta(stateSize*ii-1) = zeta(stateSize*ii-1) - QinvSM*(relHeading.deltaHeading(ii,jj)+2*pi)  ;
                    
                else
                    % diags
                    Omega(stateSize*ii-1,stateSize*ii-1) = Omega(stateSize*ii-1:stateSize*ii-1) + QinvSM;
                    Omega(stateSize*jj-1,stateSize*jj-1) = Omega(stateSize*jj-1:stateSize*jj-1) + QinvSM;
                    % offdiags
                    Omega(stateSize*ii-1,stateSize*jj-1) = Omega(stateSize*ii-1:stateSize*ii-1) - QinvSM;
                    Omega(stateSize*jj-1,stateSize*ii-1) = Omega(stateSize*jj-1:stateSize*jj-1) - QinvSM;
                    % zeta
                    zeta(stateSize*jj-1) = zeta(stateSize*jj-1) + QinvSM*relHeading.deltaHeading(ii,jj);
                    zeta(stateSize*ii-1) = zeta(stateSize*ii-1) - QinvSM*relHeading.deltaHeading(ii,jj);
                end
            end
        end
    end
end

%% Loop closure
if (~isempty(LCobj))
    fprintf('Loop closure')
    % translation
    R = Euler2RotMat(0,0,stateEstimate(5,LCobj.idx1));
    z = R(1:2,1:2)*LCobj.T_1_2(1:2);
    H = [eye(2), zeros(2,4),   -eye(2), zeros(2,4)];
    muVec = [stateEstimate(:,LCobj.idx1); stateEstimate(:,LCobj.idx2)];
    zHat = H*muVec;
    QinvX = 1e8*eye(2);
    CovAdd = H'*QinvX*H;
    zetaAdd = H'*QinvX*(z - zHat + H*muVec); % redundant, since this measurement is linear
    % Add it to the matrices
    ii = LCobj.idx1;
    jj = LCobj.idx2;
    % add info to omega
    Omega(ii*stateSize-5:ii*stateSize,ii*stateSize-5:ii*stateSize) = ... %self
        sparse(Omega(ii*stateSize-5:ii*stateSize,ii*stateSize-5:ii*stateSize)...
        + CovAdd(1:stateSize,1:stateSize));
    Omega(ii*stateSize-5:ii*stateSize,jj*stateSize-5:jj*stateSize) = ... %self
        sparse(Omega(ii*stateSize-5:ii*stateSize,jj*stateSize-5:jj*stateSize)...
        + CovAdd(1:stateSize,stateSize+1:end));
    Omega(jj*stateSize-5:jj*stateSize,ii*stateSize-5:ii*stateSize) = ... %self
        sparse(Omega(jj*stateSize-5:jj*stateSize,ii*stateSize-5:ii*stateSize)...
        + CovAdd(stateSize+1:end,1:stateSize));    
    Omega(jj*stateSize-5:jj*stateSize,jj*stateSize-5:jj*stateSize) = ... %self
        sparse(Omega(jj*stateSize-5:jj*stateSize,jj*stateSize-5:jj*stateSize)...
        + CovAdd(stateSize+1:end,stateSize+1:end));
    % Add info to zeta
    zeta((ii-1)*stateSize+1:(ii)*stateSize) = zeta((ii-1)*stateSize+1:(ii)*stateSize) +...
        + zetaAdd(1:stateSize);
    zeta((jj-1)*stateSize+1:(jj)*stateSize) = zeta((jj-1)*stateSize+1:(jj)*stateSize) +...
        + zetaAdd(stateSize+1:end);    
%     % rotation
%     z = -asin(LCobj.R_1_2(2,1)); % "measured" delta psi
%     H = [0 0 0 0 1 0 0 0 0 0 -1 0];
%     muVec = [stateEstimate(:,LCobj.idx1); stateEstimate(:,LCobj.idx2)];
%     R12hat = Euler2RotMat(0,0,stateEstimate(5,LCobj.idx1))'*Euler2RotMat(0,0,stateEstimate(5,LCobj.idx2)) ;
%     zHat = asin(R12hat(2,1));
%     QinvX = 10;
%     CovAdd = H'*QinvX*H;
%     zetaAdd = H'*QinvX*(z - zHat + mod(H*muVec,2*pi)); 
%     % add info to omega
%     Omega(ii*stateSize-5:ii*stateSize,ii*stateSize-5:ii*stateSize) = ... %self
%         sparse(Omega(ii*stateSize-5:ii*stateSize,ii*stateSize-5:ii*stateSize)...
%         + CovAdd(1:stateSize,1:stateSize));
%     Omega(ii*stateSize-5:ii*stateSize,jj*stateSize-5:jj*stateSize) = ... %self
%         sparse(Omega(ii*stateSize-5:ii*stateSize,jj*stateSize-5:jj*stateSize)...
%         + CovAdd(1:stateSize,stateSize+1:end));
%     Omega(jj*stateSize-5:jj*stateSize,ii*stateSize-5:ii*stateSize) = ... %self
%         sparse(Omega(jj*stateSize-5:jj*stateSize,ii*stateSize-5:ii*stateSize)...
%         + CovAdd(stateSize+1:end,1:stateSize));    
%     Omega(jj*stateSize-5:jj*stateSize,jj*stateSize-5:jj*stateSize) = ... %self
%         sparse(Omega(jj*stateSize-5:jj*stateSize,jj*stateSize-5:jj*stateSize)...
%         + CovAdd(stateSize+1:end,stateSize+1:end));
%     % Add info to zeta
%     zeta((ii-1)*stateSize+1:(ii)*stateSize) = zeta((ii-1)*stateSize+1:(ii)*stateSize) +...
%         + zetaAdd(1:stateSize);
%     zeta((jj-1)*stateSize+1:(jj)*stateSize) = zeta((jj-1)*stateSize+1:(jj)*stateSize) +...
%         + zetaAdd(stateSize+1:end);    
end

%% Measurements

fprintf('\nMeasurements...\n')
FeatureIndices = unique(correspondences);
FeatureIndices = FeatureIndices(2:end); % remove -17
counter = 0;
% for mapFeatureIterator = 1:length(FeatureIndices)%NmapFeatures;
%     j = FeatureIndices(mapFeatureIterator);
%     [j_s,i_s] = find(correspondences == j);
%     for idummy = 1:length(i_s)
%         ii = i_s(idummy);
%         jj = j_s(idummy);
for ii = 1:Nstates
    %fprintf('Nstates: %d\n',Nstates)
    
    B_R_V = Euler2RotMat(0,0,stateEstimate(5,ii)); 
    V_R_B = B_R_V';
    for jj = 1:sum(meas_ind(:,ii) ~= -17)
        j = correspondences(jj,ii);
            %fprintf('ii: %d, jj: %d, j: %d\n',ii, jj,j)
        measIndex = meas_ind(jj,ii); 
        mapFeatureIterator = j;
        %if(correspondences(jj,ii) == j)
            %expectedMeasurement = B_R_V'*(mapEstimate(:,j) - [stateEstimate(1:2,ii); 0] + [0 0 eps]') ;
            %expectedMeasurement = B_R_V'*(-[stateEstimate(1:2,ii); 0] + [0 0 eps]') ;
            % build covariance
            Qinv = 100*eye(3);

            xDiff = mapEstimate(1:3,j) - [stateEstimate(1:2,ii); 0];
            zHat = V_R_B*[xDiff];
            zMeas = rangeMeasurements(:,measIndex);
            cTheta = B_R_V(1,1);
            sTheta = B_R_V(2,1);
            %[zHat zMeas zDiff]
            
            %             H = [-cTheta, -sTheta, 0, 0, -xDiff(1)*sTheta+xDiff(2)*cTheta, 0, cTheta, sTheta, 0;...
            %                 sTheta, -cTheta, 0, 0, -xDiff(1)*cTheta-xDiff(2)*sTheta, 0, -sTheta, cTheta, 0;...
            %                 0,0,0,0,0,0,0,0,1];
            H = [[-V_R_B(1:2,1:2), zeros(2), [-sTheta cTheta; -cTheta -sTheta]*xDiff(1:2) , [0;0],  V_R_B(1:2,1:3)];...
                [0, 0,           0,0                   ,0                         ,0       ,0,0,  1]];
            %            H = [[-V_R_B(1:2,1:2), zeros(2,4),  V_R_B(1:2,1:3)];...
            %                     [0, 0,           0,0                   ,0                         ,0       ,0,0,  1]];
            %             H = [[-eye(2), zeros(2,4),  eye(2), [0;0] ];...
            %                     [0, 0,           0,0                   ,0                         ,0       ,0,0,  1]];
            
            CovAdd = H'*Qinv*H;
            %CovAdd(end,end) = .01;
            %zetaAdd = H'*Qinv*(H*[stateEstimate(:,ii) ; mapEstimate(:,j)]);
            zetaAdd = H'*Qinv*(zMeas - zHat + H*[stateEstimate(:,ii) ; mapEstimate(:,j)]);
            
            % Add information about position to Omega and zeta
            
            Omega(ii*stateSize-5:ii*stateSize,ii*stateSize-5:ii*stateSize) = ... %self
                sparse(Omega(ii*stateSize-5:ii*stateSize,ii*stateSize-5:ii*stateSize)...
                + CovAdd(1:stateSize,1:stateSize));
            
            Omega(ii*stateSize-5:ii*stateSize,Nstates*stateSize+mapFeatureIterator*3-2:Nstates*stateSize+(mapFeatureIterator*3)) = ... %self
                sparse(Omega(ii*stateSize-5:ii*stateSize,Nstates*stateSize+mapFeatureIterator*3-2:Nstates*stateSize+(mapFeatureIterator*3))...
                + CovAdd(1:stateSize,stateSize+1:stateSize+3));
            
            zeta((ii-1)*stateSize+1:(ii)*stateSize) = zeta((ii-1)*stateSize+1:(ii)*stateSize) +...
                + zetaAdd(1:stateSize);
            
            %          % Add information about map to Omega and zeta
            Omega(Nstates*stateSize+(mapFeatureIterator-1)*3+1:Nstates*stateSize+(mapFeatureIterator*3),Nstates*stateSize+(mapFeatureIterator-1)*3+1:Nstates*stateSize+(mapFeatureIterator*3)) = ... %self
                Omega(Nstates*stateSize+(mapFeatureIterator-1)*3+1:Nstates*stateSize+(mapFeatureIterator*3),Nstates*stateSize+(mapFeatureIterator-1)*3+1:Nstates*stateSize+(mapFeatureIterator*3))...
                + sparse(CovAdd(stateSize+1:end,stateSize+1:end));
            Omega(Nstates*stateSize+(mapFeatureIterator-1)*3+1:Nstates*stateSize+(mapFeatureIterator*3),(ii-1)*stateSize+1:(ii)*stateSize) = ...
                Omega(Nstates*stateSize+(mapFeatureIterator-1)*3+1:Nstates*stateSize+(mapFeatureIterator*3),(ii-1)*stateSize+1:(ii)*stateSize) ...
                + sparse(CovAdd(stateSize+1:stateSize+3,1:stateSize));
            
            zeta(Nstates*stateSize+(mapFeatureIterator-1)*3+1:Nstates*stateSize+(mapFeatureIterator*3)) = ...
                zeta(Nstates*stateSize+(mapFeatureIterator-1)*3+1:Nstates*stateSize+(mapFeatureIterator*3)) + ...
                + zetaAdd(stateSize+1:end);
        %end
        

    end
    if(false)%counter > 100)
        counter = 0;
        fprintf('%d seconds\n',timeStamps(ii))
    end
    counter = counter+1;
    

    
    
end
c_i_t_new = correspondences;


for ii = 1:Nstates
    % Handle dvl measurements
    % and penalize large iceberg angular velocities ?
    
%     %% omega:
    % initial idea: 2*sigma = .01, so sigma = .005. 1/.005^2 = 40000
    penalty = 100; % 1/covariance of expected measurement
    % Penalize large omega
    Homega = [0 0 0 0 0 1];
    zMeasOmega = stateEstimate(6,ii);
    zExpOmega = 0;
    Omega(ii*stateSize,ii*stateSize) = Omega(ii*stateSize,ii*stateSize) + penalty;
    zeta(ii*stateSize-5:ii*stateSize) = zeta(ii*stateSize-5:ii*stateSize) + penalty*Homega'*(zMeasOmega);
%    zeta(ii*stateSize-5:ii*stateSize) = zeta(ii*stateSize-5:ii*stateSize) + penalty*Homega'*(zMeasOmega - zExpOmega + Homega*stateEstimate(:,ii));

    %% dvl
    dvlMeas = dvlReadings(ii) ;
    Qdvl = 1e-6;
    for jDvl = 1:dvl.numBeams
        if(~isnan(dvlMeas.ranges(jDvl))) % valid return
            
            Hdvl_j = [0 0 dvl.beamsVF(1:2,jDvl)' 0 0];
            zMeasDVL = -dvlMeas.normalVelocity(jDvl);
            zHatDVL = Hdvl_j*stateEstimate(:,ii);
            CovAdd = 1/Qdvl*(Hdvl_j'*Hdvl_j);
            zetaAdd = 1/Qdvl*Hdvl_j'*(zMeasDVL - zHatDVL + Hdvl_j*stateEstimate(:,ii)); % I think this is right, because the dvl meas fxn is linear.
            Omega((ii-1)*stateSize+1:(ii)*stateSize,(ii-1)*stateSize+1:(ii)*stateSize) = ... %self
                    sparse(Omega((ii-1)*stateSize+1:(ii)*stateSize,(ii-1)*stateSize+1:(ii)*stateSize)...
                    + CovAdd);
                zeta((ii-1)*stateSize+1:(ii)*stateSize) = zeta((ii-1)*stateSize+1:(ii)*stateSize) +...
                    + zetaAdd;
            end
        end 
    
    
end


% if(haveInitialized)
%     for jjj = 1:length(remainingFeatures)
%         
%         c_i_t_new(correspondences==remainingFeatures(jjj)) = jjj;
%         if(mod(jjj,50))
%             fprintf('%d\n',jjj)
%         end
%     end
% end

end

    function deriv = dStateDt(~,stateAndInputs)
        
        head = stateAndInputs(5);
        R = Euler2RotMat(0,0,head);
        b_R_v = R(1:2,1:2);
        
        deriv = zeros(size(stateAndInputs));
        deriv(1:2) = b_R_v*stateAndInputs(3:4);
        deriv(3:4) = [stateAndInputs(7); 0] - stateAndInputs(3:4);
        deriv(5) = (stateAndInputs(8) - stateAndInputs(6));
        deriv(6) = 0;%-.00001*stateAndInputs(6);
        % zeros for everything else
    end











