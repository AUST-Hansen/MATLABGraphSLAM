function [Omega,zeta] = GraphSLAM_linearize(timeStamps,inputs,measTimestamps,rangeMeasurements,dvl,dvlReadings,correspondences,stateEstimate)

   [stateSize, Nstates] = size(stateEstimate);
   NmapObservations = size(rangeMeasurements,2);
   % state has position, velocity, and heading info, as well as berg omega
   % estimate. Map features have 3d position.
   Omega = spalloc((stateSize*Nstates+3*NmapObservations),(stateSize*Nstates+3*NmapObservations),...
           1e7);
   % anchor initial position and orientation
   Omega(1,1) = 1e13;
   Omega(2,2) = 1e13;
   Omega(3,3) = 1;
   Omega(4,4) = 1;
   Omega(5,5) = 1e13;
   Omega(6,6) = 1;
   % ininitialize zeta
   zeta = zeros(stateSize*Nstates + 3*NmapObservations,1);
   
   xhat = zeros(stateSize,1);
   G = zeros(stateSize);
   velocityDamping = .9;
   processNoise = diag([.1,.1,.1,.1,.001,.001]);

   %% Motion model
   for ii = 1:Nstates-1
       dT = (timeStamps(ii+1) - timeStamps(ii));
       heading = stateEstimate(5,ii);
       berg_R_veh = [cos(heading), -sin(heading); sin(heading), cos(heading)]; 
       xhat(1:2) = eye(2)*stateEstimate(1:2,ii) + berg_R_veh*stateEstimate(3:4,ii)*dT; % update position
       xhat(3:4) = .9*stateEstimate(3:4,ii) + .1*[inputs(1,ii); 0];
       xhat(5) = stateEstimate(5,ii) + (inputs(2,ii) - stateEstimate(6,ii))*dT;
       xhat(6) = stateEstimate(6,ii); 
       
       G(1:2,1:2) = eye(2);
       G(1:2,3:4) = berg_R_veh*dT;
       G(3:4,3:4) = velocityDamping*eye(2);
       G(5:6,5:6) = eye(2);
       G(5,6) = -dT;
       % add information to Omega and zeta
       Omega((ii-1)*stateSize+1:(ii+1)*stateSize,(ii-1)*stateSize+1:(ii+1)*stateSize) = ... %self
                  Omega((ii-1)*stateSize+1:(ii+1)*stateSize,(ii-1)*stateSize+1:(ii+1)*stateSize)...
                + [-G'; eye(stateSize)]*(processNoise\[-G, eye(stateSize)]);
            
       zeta((ii-1)*stateSize+1:(ii+1)*stateSize) = zeta((ii-1)*stateSize+1:(ii+1)*stateSize) +...
               [-G'; eye(stateSize)]*(processNoise\(xhat - G*stateEstimate(:,ii)));
           
   end

   %% Measurements
   % measurement noise
   for ii = 1:Nstates
       % iterate through all valid range measurements
       measIterator = 1;
       j = correspondences(measIterator,ii);
       while (j > 0) % failed correspondences will have negative j value
          
       mu
           
           
       % increment j    
       measIterator = measIterator+1;
       j = correspondences(measIterator,ii);
       end
   end
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   