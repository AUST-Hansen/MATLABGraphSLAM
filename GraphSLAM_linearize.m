function [Omega,zeta] = GraphSLAM_linearize(timeStamps,inputs,measTimestamps,sensor,rangeMeasurements,dvl,dvlReadings,correspondences,fullStateEstimate)

   Nstates = size(timeStamps,2);
   stateSize = 6;
   NmapObservations = size(rangeMeasurements,2);
   stateEstimate = reshape(fullStateEstimate(1:stateSize*Nstates),stateSize,[]);
   mapEstimate = reshape(fullStateEstimate(stateSize*Nstates+1:end),3,[]);
   % state has position, velocity, and heading info, as well as berg omega
   % estimate. Map features have 3d position.
   Omega = spalloc((stateSize*Nstates+3*NmapObservations),(stateSize*Nstates+3*NmapObservations),10e7);
   % anchor initial position and orientation
   Omega(1,1) = 1e13;
   Omega(2,2) = 1e13;
   Omega(3,3) = 1;
   Omega(4,4) = 1e13;
   Omega(5,5) = 1e13;
   Omega(6,6) = 1e13;
   % ininitialize zeta
   zeta = zeros(stateSize*Nstates + 3*NmapObservations,1);
   zeta(1:6) = [0 0 inputs(1,1) 0 inputs(2,1) 0]';
   xhat = zeros(stateSize,1);
   G = zeros(stateSize);
   velocityDamping = .1;
   processNoise = diag([.1,.1,.1,.1,.001,.001]);

   %% Motion model
   for ii = 1:Nstates-1
       dT = (timeStamps(ii+1) - timeStamps(ii));
       heading = stateEstimate(5,ii);
       berg_R_veh = [cos(heading), -sin(heading); sin(heading), cos(heading)]; 
       xhat(1:2) = eye(2)*stateEstimate(1:2,ii) + berg_R_veh*stateEstimate(3:4,ii)*dT; % update position
       xhat(3:4) = .0*stateEstimate(3:4,ii) + 1*[inputs(1,ii); 0];
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
   counter = 0;
   for ii = 1:Nstates
       % iterate through all valid range measurements
       measIterator = 1;
       j = correspondences(measIterator,ii);
       while (j > 0) % failed correspondences will have negative j value
          delta = mapEstimate(:,j) - [stateEstimate(1:2,ii); 0] + [0 0 eps]'; % switch sign!!!!!!!!!!!!!!!!!!!!
          q = delta'*delta;
          sqrtQ = sqrt(q)+eps;
          Qrange = diag([(sqrtQ*sensor.beamSigmaPercentageOfRange)^2 .001 .001]);
          Qinv = inv(Qrange);
          bearingHat = atan2(delta(2),delta(1));
          if bearingHat < 0
              bearingHat = bearingHat + 2*pi;
          end
          currentHeading = stateEstimate(5,ii);
          while(currentHeading > 2*pi)
              currentHeading = currentHeading - 2*pi;
          end
          while(currentHeading < 0)
              currentHeading = currentHeading + 2*pi;
          end
          expBearing = bearingHat - currentHeading;
          while (expBearing >= 2*pi)
              expBearing = expBearing - 2*pi;
          end
          while (expBearing <= 0)
              expBearing = expBearing + 2*pi;
          end
          zHat = [sqrtQ; expBearing; asin(delta(3)/sqrtQ)];
          zMeas = rangeMeasurements(:,j);
          zDiff = zMeas - zHat;
          %[zHat zMeas zDiff]
          if(zDiff(2) < -pi)
              zDiff(2) = zDiff(2) + 2*pi;
          elseif (zDiff(2) > pi)
              zDiff(2) = zDiff(2) - 2*pi;
          end
           H = 1/q*[-delta(1) -delta(2) 0 0 0 0 delta(1) delta(2) delta(3);...
                    delta(2) -delta(1) 0 0 -q 0 -delta(2) delta(1)  0;...
                    0 0 0 0 0 0 0 0 1/sqrt(1-(delta(3)/sqrtQ)^2)*(2*delta(3)) ];
%           Hx = 1/q*[-sqrtQ*delta(1) -sqrtQ*delta(2) 0 0 0 0;...
%                    delta(2) -delta(1) 0 0 -q 0 ;...
%                    0 0 0 0 0 0];
%           Hm = 1/q*[sqrtQ*delta(1) sqrtQ*delta(2) sqrtQ*delta(3);...
%                    -delta(2) delta(1)  0; ...
%                    0 0 1/sqrt(1-(delta(3)/sqrtQ)^2)*(2*delta(3)) ];
          CovAdd = H'*Qinv*H;
          CovAdd(end,end) = .01;
          zetaAdd = H'*Qinv*(zDiff  + H*[stateEstimate(:,ii) ; mapEstimate(:,j)]);
          
         % Add information about position to Omega and zeta

          Omega((ii-1)*stateSize+1:(ii)*stateSize,(ii-1)*stateSize+1:(ii)*stateSize) = ... %self
                   sparse(Omega((ii-1)*stateSize+1:(ii)*stateSize,(ii-1)*stateSize+1:(ii)*stateSize)...
                 + CovAdd(1:stateSize,1:stateSize));
          Omega((ii-1)*stateSize+1:(ii)*stateSize,Nstates*stateSize+(j-1)*3+1:Nstates*stateSize+(j*3)) = ... %self
                   sparse(Omega((ii-1)*stateSize+1:(ii)*stateSize,Nstates*stateSize+(j-1)*3+1:Nstates*stateSize+(j*3))...
                 + CovAdd(1:stateSize,stateSize+1:stateSize+3));
          zeta((ii-1)*stateSize+1:(ii)*stateSize) = zeta((ii-1)*stateSize+1:(ii)*stateSize) +...
                 + zetaAdd(1:stateSize);
             
%          % Add information about map to Omega and zeta
         Omega(Nstates*stateSize+(j-1)*3+1:Nstates*stateSize+(j*3),Nstates*stateSize+(j-1)*3+1:Nstates*stateSize+(j*3)) = ... %self
                Omega(Nstates*stateSize+(j-1)*3+1:Nstates*stateSize+(j*3),Nstates*stateSize+(j-1)*3+1:Nstates*stateSize+(j*3))...
                + sparse(CovAdd(stateSize+1:end,stateSize+1:end));
         Omega(Nstates*stateSize+(j-1)*3+1:Nstates*stateSize+(j*3),(ii-1)*stateSize+1:(ii)*stateSize) = ...   
                Omega(Nstates*stateSize+(j-1)*3+1:Nstates*stateSize+(j*3),(ii-1)*stateSize+1:(ii)*stateSize) ...
                + sparse(CovAdd(stateSize+1:stateSize+3,1:stateSize));
        
          zeta(Nstates*stateSize+(j-1)*3+1:Nstates*stateSize+(j*3)) = ...
              zeta(Nstates*stateSize+(j-1)*3+1:Nstates*stateSize+(j*3)) + ...
              + zetaAdd(stateSize+1:end);
         % increment j    
         measIterator = measIterator+1;
         j = correspondences(measIterator,ii);
       end
       % Handle dvl measurements
       dvlMeas = dvlReadings(ii) ;      
       Qdvl = .01;
       for jDvl = 1:dvl.numBeams
          if(~isnan(dvlMeas.ranges(jDvl))) % valid return 
           
              Hdvl_j = [0 0 dvl.beamsVF(1:2,jDvl)' 0 0];
              zMeasDVL = -dvlMeas.normalVelocity(jDvl);
              zHatDVL = Hdvl_j*stateEstimate(:,ii);
              CovAdd = 1/Qdvl*Hdvl_j'*Hdvl_j;
              zetaAdd = 1/Qdvl*Hdvl_j'*(zMeasDVL - zHatDVL + Hdvl_j*stateEstimate(:,ii)); % I think this is right, because the dvl meas fxn is linear.
              Omega((ii-1)*stateSize+1:(ii)*stateSize,(ii-1)*stateSize+1:(ii)*stateSize) = ... %self
                   sparse(Omega((ii-1)*stateSize+1:(ii)*stateSize,(ii-1)*stateSize+1:(ii)*stateSize)...
                 + CovAdd);
             zeta((ii-1)*stateSize+1:(ii)*stateSize) = zeta((ii-1)*stateSize+1:(ii)*stateSize) +...
                 + zetaAdd;
          end
       end
       
       if(counter > 100)
           counter = 0;
           fprintf('%d seconds\n',timeStamps(ii))
       end
       counter = counter+1;
   end
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   