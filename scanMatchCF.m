close all;


[bodyvel, lock] = processDVL(dvl,dvlData,true,false);


% state to be estimated = [x y heading omegaberg]
measurements = zeros(1,size(omega_vehicle,2));
errors = measurements;
resids = measurements;
wignx = measurements;
state = zeros(4,size(omega_vehicle,2));
step = poseSkip;
N = size(omega_vehicle,2);
for ii = 2:1:size(omega_vehicle,2)
    ii
    % integrated term
    I_R_B = Euler2RotMat(0,0,state(3,ii-1));
    state(1:2,ii) = state(1:2,ii-1) + I_R_B(1:2,1:2)*bodyvel(1:2,ii-1)*(timeSteps(ii)-timeSteps(ii-1));
    state(3,ii) = state(3,ii-1) + (omega_vehicle(ii-1) - state(4,ii-1))*(timeSteps(ii)-timeSteps(ii-1));
    %state(3,ii) = state(3,ii-1) + (omega_vehicle(ii-1) - (euler_berg_t(3,ii)-euler_berg_t(3,ii-1))*.5) * (timeSteps(ii)-timeSteps(ii-1));
    state(4,ii) = state(4,ii-1);
    measurements(ii) = measurements(ii-1) + omega_vehicle(ii-1)*(timeSteps(ii)-timeSteps(ii-1));
    errors(ii) = errors(ii-1);
    wignx(ii) = wignx(ii-1);
    % Scanmatching
    idx1 = ii-step;
    idx2 = ii;
    try
        pointCloud1 = rangeMeasurements(:,meas_ind(meas_ind(:,idx1)~=-17,idx1));
        pointCloud2 = rangeMeasurements(:,meas_ind(meas_ind(:,idx2)~=-17,idx2));
    catch
        pointCloud1 = []
    end
    
    if(~isempty(pointCloud1))
        [TR, TT, ER, knnOut] = robustScanMatch(pointCloud1,pointCloud2);
        
        if (~isempty(TR))
            resids(ii) = ER;
            yawMeas = asin(-TR(1,2)); %!!!!!!!!!!!!!!!!!!! SIGN?
            180/pi*yawMeas;
            
            yawScanMatch = state(3,idx1) + yawMeas;
            
            
            measurements(ii) = yawScanMatch;
            integratedOmega = sum(omega_vehicle(idx1:idx2) - state(4,idx1:idx2))*.5;
            wignx(ii) = yawMeas/(timeSteps(idx2)-timeSteps(idx1));
            %errors(ii) = .7*errors(ii) + .3*(yawMeas - integratedOmega)/(timeSteps(idx2)-timeSteps(idx1));
            errors(ii) = +(yawMeas - integratedOmega);
            % valid measurement: update
            alfa1 = .95;
            alfa2 = .01;
            state(4,ii) = state(4,ii) + alfa2*max(min(errors(ii),1e-3),-1e-3);
            state(3,ii) = alfa1*state(3,ii) + (1-alfa1)*measurements(ii);
            
            
        end
        
    end
    
    %Input = input('Press enter to continue\n','s');
    
end
%%
figure;
plot(-state(2,:),state(1,:));
axis equal
hold on
plot(xVehBergframe(1,:)- xVehBergframe(1,1),xVehBergframe(2,:)- xVehBergframe(2,1),'g');

figure; 
subplot(3,1,1)
plot(state(3,:)); hold on; 
plot(measurements,'g'); 
plot(euler_obs_t(3,1:N) - euler_berg_t(3,1:N) - pi/2,'r.')
plot(euler_obs_t(3,1:N) -pi/2,'r')
plot(10*euler_berg_t(3,1:N),'k')
legend('estimated heading','imagenex measurements','bergframe heading','inertial heading','10*iceberg heading')

subplot(3,1,2)
omegaBerg = diff(euler_berg_t(3,1:N))/.5;
avgOmega = mean(omegaBerg)*ones(size(omegaBerg));
plot(omegaBerg)
hold on;
%plot(avgOmega,'r')
plot(state(4,:),'g')
DEL = diag(-ones(1,N)) + diag(ones(1,N-step),step);
diffs = DEL*measurements';
plot(smooth(omega_vehicle-wignx,30),'k')
integ = 0*(DEL*measurements');
for ii = 2:N
    integ(ii) = integ(ii-1) + diffs(ii);
end
legend('wBerg actual','estimated omega','wignx')
subplot(3,1,3)
figure;
plot(omega_vehicle - state(4,:));
hold on
plot(omega_vehicle - [0 omegaBerg],'r');
plot(omega_vehicle,'g') 
plot(smooth(wignx,30),'c')
plot(.01*resids,'k')
legend('berg w veh est','berg w veh act','inertial omega','wignx','resids')