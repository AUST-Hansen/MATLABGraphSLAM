%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    GraphSLAMicebergScript.m
%
%    Runs GraphSLAM as described in pp. 347-365 of Probabilistic Robotics
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
pause on
%% Constant parameters
% File to read
filename = 'RandomTrajectories/WallFollowingTraj2.mat';
%filename = 'RandomTrajectories/WideAngleWallFollowingTraj1.mat';
% Number of points to use
% startIdx and endIdx will be used for all data. skip will only be applied to
% range data, to allow a sparser map of observed features.
startIdx = 1; % for now, this must be 1
skip = 1;     % for now, this must be 1 also
rangeSkip = 1;
poseSkip = 10;
endIdx = 800;
stateSize = 6;
% Imagenex matching parameters
ignx_sparsity = 3;
scanmatch_threshold = 20; % only look for matches within this many meters of pose difference
scanmatch_RMStolerance = 2.5;
scanmatch_probThreshold = 1; 
% Use imagenex, multibeam and DVL ranges in observations?
useDVLRanges = false;
useMultibeamRanges = false;
useImagenexData = true;
% give loop closure a few iterations to do its thing
lcAllowance = 2;
MAX_ITER = 40;
%% Read in Data
fprintf('Reading in Data...\n')
addpath ..
load(filename)
%% Gather up truth for comparison
fprintf('Processing "Truth" data...\n')
fprintf('Imagenex...\n')
[trueIgnxCloud, xVehBergframe] = reconstructTrueCloud(x_obs_t_true(:,startIdx:skip:endIdx),euler_obs_t(:,startIdx:skip:endIdx),x_berg_cm_t(:,startIdx:skip:endIdx),euler_berg_t(:,startIdx:skip:endIdx),imagenex,imagenexData(:,startIdx:skip:endIdx));
trueIgnxCloud = trueIgnxCloud - repmat(x_obs_t_true(:,1),1,size(trueIgnxCloud,2));
fprintf('Reson...\n')
[trueBergShape,~] = reconstructTrueCloud(x_obs_t_true(:,startIdx:skip:endIdx),euler_obs_t(:,startIdx:skip:endIdx),x_berg_cm_t(:,startIdx:skip:endIdx),euler_berg_t(:,startIdx:skip:endIdx),sensor,truthrangeData(:,startIdx:skip:endIdx));
%% Clean data of bogus ranges %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if range == -17 or MAX_Range, remove it.
%--------------------------------------------------------------------------
% MeasurementTimeStamps allows subsampling of range data, meaning we can
% integrate high rate pose data but not deal with as many map features if
% we don't want to.
fprintf('Cleaning and formatting data...\n')
[measurementTimestamps, rangeMeasurements, c_i_t, meas_ind] = cleanMeasurements(timeSteps(1:skip:endIdx),sensor,rangeData(:,startIdx:skip:endIdx),useMultibeamRanges,imagenex,imagenexData(:,startIdx:skip:endIdx),useImagenexData,dvl,dvlData(startIdx:skip:endIdx),useDVLRanges,rangeSkip,poseSkip);
%% try to identify loop closure
c_i_t = lookForLoopClosure(c_i_t,rangeMeasurements,measurementTimestamps,scanmatch_probThreshold);
%--------------------------------------------------------------------------
%% Form input vectors  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract angular velocity of vehicle in inertial as a control input
% 2xendIdxtimeSteps matrix, where first row is vx_commanded and 2nd row is omega
% of vehicle in inertial
fprintf('Building inputs...\n')
v_vehicle_inertial = sqrt(sum((diff(x_obs_t(:,startIdx:endIdx)')/.5).^2,2))';    % Clean this up later 
v_commanded = [v_vehicle_inertial, v_vehicle_inertial(end)];
omega_vehicle = diff(euler_obs_t(3,startIdx:endIdx))./diff(timeSteps(startIdx:endIdx));
omega_vehicle = [omega_vehicle, omega_vehicle(end)]; % not really needed, but makes omega same length as timeSteps
inputs = [v_commanded; omega_vehicle];
%% Run GraphSLAM
fprintf('Begin GraphSLAM...\n')
fprintf('Initializing...\n')
initialPsiDotGuess = .0;%035;
initialStateEstimate = GraphSLAM_initialize(timeSteps(startIdx:endIdx),inputs,initialPsiDotGuess);
%initialStateEstimate = GraphSLAM_initializeWithTruth(timeSteps(startIdx:endIdx),euler_obs_t,euler_berg_t,x_obs_t_true,x_berg_cm_t)
fprintf('Initializing Map Features...\n')
initialMapEstimate = GraphSLAM_initializeMap(initialStateEstimate,rangeMeasurements,c_i_t);
initialFullStateEstimate = [reshape(initialStateEstimate,[],1); reshape(initialMapEstimate,[],1)];
%% Calculate inzformation form of full posterior
fprintf('Linearizing...\n')
[Omega,zeta] = GraphSLAM_linearize(timeSteps(startIdx:endIdx),inputs,measurementTimestamps,imagenex,rangeMeasurements,dvl,dvlData,c_i_t,meas_ind,initialFullStateEstimate,false);
%% initial testing stuff
%stateHist = reshape(Omega(1:endIdx*6,1:endIdx*6)\zeta(1:endIdx*6),6,[]);
figure(1);
plot(timeSteps,euler_berg_t(3,:));
title('iceberg heading')
state = Omega\zeta;
stateHist = reshape(state(1:6*(endIdx-startIdx)),6,[]);
mapEsts = reshape(state(6*endIdx+1:end),3,[]);
figure(2)
plot(-initialStateEstimate(2,:),initialStateEstimate(1,:))
hold on
axis equal
plot(-stateHist(2,:),stateHist(1,:),'r')
plot(xVehBergframe(1,:) - xVehBergframe(1,1)  ,xVehBergframe(2,:) - xVehBergframe(2,1),'g')
scatter(trueIgnxCloud(1,:),trueIgnxCloud(2,:),ones(1,size(trueIgnxCloud,2)),'g')
scatter(-initialMapEstimate(2,:),initialMapEstimate(1,:),ones(1,size(initialMapEstimate,2)),'b')
scatter(-mapEsts(2,:),mapEsts(1,:),ones(1,size(mapEsts,2)),'r')
%legend('initial estimate','after being pushed through information form','truth')
title('estimated path')
drawnow;
%% visualizing icp stuff
% first two scans
% figure;
% scatter(initialMapEstimate(1,c_i_t(c_i_t(:,1)~=-17,1)),initialMapEstimate(2,c_i_t(c_i_t(:,1)~=-17,1)),'r')
% hold on; axis equal;
% scatter(initialMapEstimate(1,c_i_t(c_i_t(:,20)~=-17,20)),initialMapEstimate(2,c_i_t(c_i_t(:,20)~=-17,20)),'b')
% title('scans 1 and 20')
% figure;
% scatter(initialMapEstimate(1,c_i_t(c_i_t(:,end-20)~=-17,end-20)),initialMapEstimate(2,c_i_t(c_i_t(:,end-20)~=-17,end-20)),'r')
% hold on; axis equal;
% scatter(initialMapEstimate(1,c_i_t(c_i_t(:,end)~=-17,end)),initialMapEstimate(2,c_i_t(c_i_t(:,end)~=-17,end)),'b')
% title('scans end-20 and end')
%%
% qICP = initialMapEstimate(:,c_i_t(c_i_t(:,end-20)~=-17,end-20));
% pICP = initialMapEstimate(:,c_i_t(c_i_t(:,end)~=-17,end));
% [R,T,ERR] = icp(qICP,pICP,'twoDee',false);
% figure;
% pNew = R*pICP + repmat(T,1,size(pICP,2));
% scatter(qICP(1,:),qICP(2,:),'r')
% hold on; axis equal;
% scatter(pNew(1,:),pNew(2,:),'b')
%% Reduce graph by marginalizing 
fprintf('Reducing...\n')
[OmegaReduced,zetaReduced] = GraphSLAM_reduce(timeSteps(startIdx:endIdx),stateSize,Omega,zeta,c_i_t);
% state2 = OmegaReduced\zetaReduced;
% stateHist2 = reshape(state2,6,[]);
% figure(1)
% plot(-stateHist2(2,:),stateHist2(1,:),'.k')
% hold on
% axis equal
% title('reduced solution')
%% Initial solve
fprintf('Solving...\n')
[mu, Sigma] = GraphSLAM_solve(OmegaReduced,zetaReduced,Omega,zeta,c_i_t);

figure(12)
plot(xVehBergframe(1,:) - xVehBergframe(1,1)  ,xVehBergframe(2,:) - xVehBergframe(2,1),'g')
hold on; axis equal
plot(-stateHist(2,:),stateHist(1,:),'k')
scatter(-mapEsts(2,:),mapEsts(1,:),ones(1,size(mapEsts,2)),'k')
stateHist = reshape(mu(1:6*endIdx),6,[]);
mapEsts = reshape(mu(6*endIdx+1:end),3,[]);

% Timeout counter for main loop
%%
itimeout = 0;
DEBUGflag = false;
c_i_t_last = c_i_t;
figure(10);spy(Omega)


while (itimeout < MAX_ITER) 
    
    % Test correspondences for imagenex scan matching.
    c_i_t_last = c_i_t;
    mu_last = mu;
    if itimeout>lcAllowance
        %c_i_t = meas_ind; % TODO: only reset once
        [c_i_t] = GraphSLAMcorrespondenceViaScanMatch(mu,c_i_t,scanmatch_threshold,scanmatch_probThreshold,scanmatch_RMStolerance,meas_ind,rangeMeasurements);
    end
        % linearize
    fprintf('Linearizing...\n')
    clear Omega
    clear Sigma
    [Omega,zeta,c_i_t] = GraphSLAM_linearize(timeSteps(startIdx:endIdx),inputs,measurementTimestamps,imagenex,rangeMeasurements,dvl,dvlData,c_i_t,meas_ind,mu,true);
    % reduce
%%
    fprintf('Reducing...\n')
    [OmegaReduced,zetaReduced] = GraphSLAM_reduce(timeSteps(startIdx:endIdx),stateSize,Omega,zeta,c_i_t);

    % solve
    fprintf('Solving...\n')
    [mu, Sigma] = GraphSLAM_solve(OmegaReduced,zetaReduced,Omega,zeta,c_i_t);
    %sanityCheck = Omega\zeta;
    %mu = .1*mu_last+.9*mu;
   itimeout=itimeout+1 ;
   fprintf('%d iterations completed...\n',itimeout)
   
   figure(10);subplot(2,1,1); spy(Omega); subplot(2,1,2); spy(OmegaReduced);
   %keyboard;
   
   if (DEBUGflag)
       
      fig3=figure(3);
      set(fig3,'Position',[0 0 700 500]);
        
        hold on
        axis equal
        stateHist = reshape(mu(1:6*(endIdx-startIdx)),6,[]);
        mapEsts = reshape(mu(6*endIdx+1:end),3,[]);
        %sanityStates = reshape(sanityCheck(1:6*(endIdx-startIdx)),6,[]);
        %sanityMaps = reshape(sanityCheck(6*endIdx+1:end),3,[]);
        colorz = 'rbkcy';
        plot(-stateHist(2,:),stateHist(1,:),colorz(mod(itimeout,5)+1))
        %plot(-sanityStates(2,:),sanityStates(1,:),'c')
        if(itimeout ==1)
          plot(xVehBergframe(1,:) - xVehBergframe(1,1)  ,xVehBergframe(2,:) - xVehBergframe(2,1),'g')
          scatter(trueIgnxCloud(1,:),trueIgnxCloud(2,:),ones(1,size(trueIgnxCloud,2)),'g')
        end
        scatter(-mapEsts(2,:),mapEsts(1,:),ones(1,size(mapEsts,2)),'k')
          for qq = 1:5:endIdx-1
              B_R_Vi = Euler2RotMat(0,0,stateHist(5,qq));
              velocity = B_R_Vi*[stateHist(3:4,qq);0];
              scatter(-stateHist(2,qq),stateHist(1,qq),colorz(mod(qq,5)+1))
              quiver(-stateHist(2,qq),stateHist(1,qq),-velocity(2),velocity(1))
              plotFeatures = find(c_i_t_last(:,qq)~=-17);
%               if(~isempty(plotFeatures))
%                 scatter(-mapEsts(2,c_i_t_last(plotFeatures,qq)),mapEsts(1,c_i_t_last(plotFeatures,qq)), ones(size(mapEsts(1,c_i_t_last(plotFeatures,qq)))))
%                 drawnow()
%               end
         end
         figure(4); plot(stateHist(3:end,:)')
         hold on;
         plot(euler_obs_t(3,1:endIdx) - euler_berg_t(3,1:endIdx) - pi/2)
   end
   if(true)
        colorz = 'rbkcy';
        stateHist = reshape(mu(1:6*(endIdx-startIdx)),6,[]);
        figure(5); plot(stateHist(end,:)',colorz(mod(itimeout,5)+1))
        hold on; plot(diff(euler_berg_t(3,1:endIdx))./diff(timeSteps(1:endIdx)),'g')
        title('iceberg angular velocity')
        legend('estimated','actual');
        figure(6); spy(Omega)
       drawnow()
       
   end
   %%
   figure(7); spy(c_i_t-c_i_t_last);
    if ( sum(sum(c_i_t - c_i_t_last)) == 0 && itimeout >lcAllowance+1)
        noNewFeatures = noNewFeatures+1;
    else
        noNewFeatures = 0;
    end
    if ( noNewFeatures >= 4|| itimeout >=MAX_ITER) 
    %if (itimeout >=12)
        figure(3)
%         for zz = 1:size(c_i_t,1)*size(c_i_t,2)
%             if (c_i_t(zz) ~= -17)
%                 scatter(-mapEsts(2,c_i_t(zz)),mapEsts(1,c_i_t(zz)), 2)   
%                 hold on;
%             end
%         end
        fprintf('No new correspondences detected! Exiting loop...\n');
        break;
    end
   %keyboard
   %c_i_t = c_i_t_last;
end

%% Print a whole buncha crap
      fig3=figure(3);
      set(fig3,'Position',[0 0 700 500]);
        
        hold on
        axis equal
        stateHist = reshape(mu(1:6*(endIdx-startIdx)),6,[]);
        mapEsts = reshape(mu(6*endIdx+1:end),3,[]);
        %sanityStates = reshape(sanityCheck(1:6*(endIdx-startIdx)),6,[]);
        %sanityMaps = reshape(sanityCheck(6*endIdx+1:end),3,[]);
        colorz = 'rbkcy';
        plot(-stateHist(2,:),stateHist(1,:),colorz(mod(itimeout,5)+1))
        %plot(-sanityStates(2,:),sanityStates(1,:),'c')
        if(itimeout ==1)
          plot(xVehBergframe(1,:) - xVehBergframe(1,1)  ,xVehBergframe(2,:) - xVehBergframe(2,1),'g')
          scatter(trueIgnxCloud(1,:),trueIgnxCloud(2,:),ones(1,size(trueIgnxCloud,2)),'g')
        end
        scatter(-mapEsts(2,:),mapEsts(1,:),ones(1,size(mapEsts,2)),'k')
          for qq = 1:5:endIdx-1
              B_R_Vi = Euler2RotMat(0,0,stateHist(5,qq));
              velocity = B_R_Vi*[stateHist(3:4,qq);0];
              scatter(-stateHist(2,qq),stateHist(1,qq),colorz(mod(qq,5)+1))
              quiver(-stateHist(2,qq),stateHist(1,qq),-velocity(2),velocity(1))
              plotFeatures = find(c_i_t_last(:,qq)~=-17);
%               if(~isempty(plotFeatures))
%                 scatter(-mapEsts(2,c_i_t_last(plotFeatures,qq)),mapEsts(1,c_i_t_last(plotFeatures,qq)), ones(size(mapEsts(1,c_i_t_last(plotFeatures,qq)))))
%                 drawnow()
%               end
         end
         figure(4); plot(stateHist(3:end,:)')
         hold on;
         plot(euler_obs_t(3,1:endIdx) - euler_berg_t(3,1:endIdx) - pi/2,'g')
         [vels,lock] = processDVL(dvl,dvlData(1:endIdx),true,false);
         plot(vels(1,:),'k')
         plot(vels(2,:),'k')
         plot(lock,'r')
         legend('v_x est','v_y est','heading est','\omega_{berg} est','true heading','DVL_x','DVL_y','DVL lock')
         title('Motion estimates vs. truth')
        figure(5); plot(stateHist(end,:)',colorz(mod(itimeout,5)+1))
        hold on; plot(diff(euler_berg_t(3,1:endIdx))./diff(timeSteps(1:endIdx)),'g')
        title('iceberg angular velocity')
        legend('estimated','actual');
        figure(6); spy(Omega)
       drawnow()

%stateHist = reshape(Omega(1:endIdx*6,1:endIdx*6)\zeta(1:endIdx*6),6,[]);
%% 
figure(3)
%spy(Omega(1:10000,1:10000)); drawnow;
%mu = Omega\zeta;
%%
state = mu;
save('OmegaAndZeta.mat','Omega','zeta');
%clear Omega;
%clear zeta;
%state = Omega\zeta;
stateHist = reshape(state(1:6*endIdx),6,[]);
mapEsts = reshape(state(6*endIdx+1:end),3,[]);
figure()
plot(xVehBergframe(1,:) - xVehBergframe(1,1)  ,xVehBergframe(2,:) - xVehBergframe(2,1),'g')
hold on; axis equal
scatter(trueIgnxCloud(1,:),trueIgnxCloud(2,:),ones(1,size(trueIgnxCloud,2)),'g')
plot(-stateHist(2,:),stateHist(1,:),'k')
scatter(-mapEsts(2,:),mapEsts(1,:),ones(1,size(mapEsts,2)),'k')
legend('initial estimate','after being pushed through information form','truth','after one iteration of correspondence')
title('estimated path')

%%

[R,T,ERR] = icp(trueIgnxCloud(:,1:10:end),[-mapEsts(2,:);mapEsts(1,:);mapEsts(3,:)],'WorstRejection',.1);
rotatedMap = R*mapEsts + repmat(T,1,size(mapEsts,2));
rotatedTrack = R(1:2,1:2)*stateHist(1:2,:) + repmat(T(1:2),1,size(stateHist,2));
%%
figure(15); 

trueIgnxCloudBeg = trueIgnxCloud(:,meas_ind(meas_ind(:,1)~=-17));

trueIgnxCloudEnd = trueIgnxCloud(:,meas_ind(meas_ind(:,741)~=-17));

plot(xVehBergframe(1,:) - xVehBergframe(1,1)  ,xVehBergframe(2,:) - xVehBergframe(2,1),'g')
hold on; axis equal
scatter(trueIgnxCloud(1,:),trueIgnxCloud(2,:),ones(1,size(trueIgnxCloud,2)),'g')
scatter(trueIgnxCloudBeg(1,:),trueIgnxCloudBeg(2,:),'ro')
scatter(trueIgnxCloudEnd(1,:),trueIgnxCloudEnd(2,:),'b+')

plot(-rotatedTrack(2,:),rotatedTrack(1,:),'k')
scatter(-rotatedMap(2,:),rotatedMap(1,:),ones(1,size(rotatedMap,2)),'k')
legend('initial estimate','after being pushed through information form','truth','after one iteration of correspondence')
title('estimated path')
