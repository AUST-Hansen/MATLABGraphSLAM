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
rangeSkip = 3;
poseSkip = 5;
endIdx = 200;
stateSize = 6;
% Imagenex matching parameters
ignx_sparsity = 3;
scanmatch_threshold = 15; % only look for matches within this many meters of pose difference
scanmatch_RMStolerance = 2;
scanmatch_probThreshold = 3; 
% Use imagenex, multibeam and DVL ranges in observations?
useDVLRanges = false;
useMultibeamRanges = false;
useImagenexData = true;

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
[measurementTimestamps, rangeMeasurements, c_i_t] = cleanMeasurements(timeSteps(1:skip:endIdx),sensor,rangeData(:,startIdx:skip:endIdx),useMultibeamRanges,imagenex,imagenexData(:,startIdx:skip:endIdx),useImagenexData,dvl,dvlData(startIdx:skip:endIdx),useDVLRanges,rangeSkip,poseSkip);
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
fprintf('Initializing Map Features...\n')
initialMapEstimate = GraphSLAM_initializeMap(initialStateEstimate,rangeMeasurements,c_i_t);
initialFullStateEstimate = [reshape(initialStateEstimate,[],1); reshape(initialMapEstimate,[],1)];
%% Calculate inzformation form of full posterior
fprintf('Linearizing...\n')
[Omega,zeta] = GraphSLAM_linearize(timeSteps(startIdx:endIdx),inputs,measurementTimestamps,imagenex,rangeMeasurements,dvl,dvlData,c_i_t,initialFullStateEstimate,false);
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
%[OmegaReduced,zetaReduced] = GraphSLAM_reduce(timeSteps(startIdx:endIdx),stateSize,Omega,zeta);
% state2 = OmegaReduced\zetaReduced;
% stateHist2 = reshape(state2,6,[]);
% figure(1)
% plot(-stateHist2(2,:),stateHist2(1,:),'.k')
% hold on
% axis equal
% title('reduced solution')
%% Initial solve
fprintf('Solving...\n')
    OmegaReduced = 0; zetaReduced = 0;
[mu, Sigma] = GraphSLAM_solve(OmegaReduced,zetaReduced,Omega,zeta,c_i_t);
% Timeout counter for main loop
%%
itimeout = 0;
MAX_ITER = 20;
DEBUGflag = true;
c_i_t_last = c_i_t;
while (itimeout < MAX_ITER) 
    
    % Test correspondences for imagenex scan matching.
    c_i_t_last = c_i_t;
    [c_i_t] = GraphSLAMcorrespondenceViaScanMatch(mu,c_i_t,scanmatch_threshold,scanmatch_probThreshold,scanmatch_RMStolerance);
    % linearize
    fprintf('Linearizing...\n')
    clear Omega
    clear Sigma
    [Omega,zeta,c_i_t] = GraphSLAM_linearize(timeSteps(startIdx:endIdx),inputs,measurementTimestamps,imagenex,rangeMeasurements,dvl,dvlData,c_i_t,mu,true);
    % reduce
%%
    fprintf('Reducing...\n')
%    [OmegaReduced,zetaReduced] = GraphSLAM_reduce(timeSteps(startIdx:endIdx),stateSize,Omega,zeta);

    % solve
    fprintf('Solving...\n')
    [mu, Sigma] = GraphSLAM_solve(OmegaReduced,zetaReduced,Omega,zeta,c_i_t);
    
    if (false)%sum(sum(c_i_t - c_i_t_last)) == 0)
        fprintf('No new correspondences detected! Exiting loop...\n');
        break;
    end
   itimeout=itimeout+1 ;
   fprintf('%d iterations completed...\n',itimeout)
   
   if (DEBUGflag)
       
      fig3=figure(3);
      set(fig3,'Position',[0 0 700 500]);
        plot(xVehBergframe(1,:) - xVehBergframe(1,1)  ,xVehBergframe(2,:) - xVehBergframe(2,1),'g')
        hold on
        axis equal
        stateHist = reshape(mu(1:6*(endIdx-startIdx)),6,[]);
        mapEsts = reshape(mu(6*endIdx+1:end),3,[]);
        colorz = 'rbkcy';
        plot(-stateHist(2,:),stateHist(1,:),colorz(mod(itimeout,5)+1))
        inx = unique(c_i_t);
        for qq = 1:endIdx
            plotFeatures = find(c_i_t_last(:,qq)~=-17);
           scatter(-mapEsts(2,c_i_t_last(plotFeatures,qq)),mapEsts(1,c_i_t_last(plotFeatures,qq)), 3*ones(size(mapEsts(1,c_i_t_last(plotFeatures,qq)))))
           drawnow()
        end
         figure(4); plot(stateHist(3:end,:)')
        figure(5); plot(stateHist(end,:)',colorz(mod(itimeout,5)+1))
        hold on; plot(diff(euler_berg_t(3,1:endIdx))./diff(timeSteps(1:endIdx)),'g')
        title('iceberg angular velocity')
        legend('estimated','actual');
        figure(6); spy(Omega)
       drawnow()
       
   end
   keyboard
   c_i_t = c_i_t_last;
end

%stateHist = reshape(Omega(1:endIdx*6,1:endIdx*6)\zeta(1:endIdx*6),6,[]);
%% 
figure(3)
%spy(Omega(1:10000,1:10000)); drawnow;
mu = Omega\zeta;
%%
state = mu;
save('OmegaAndZeta.mat','Omega','zeta');
clear Omega;
clear zeta;
stateHist = reshape(state(1:6*endIdx),6,[]);
mapEsts = reshape(state(6*endIdx+1:end),3,[]);
figure(2)
plot(-stateHist(2,:),stateHist(1,:),'k')
%scatter(-mapEsts(2,:),mapEsts(1,:),ones(1,size(mapEsts,2)),'k')
legend('initial estimate','after being pushed through information form','truth','after one iteration of correspondence')
title('estimated path')
