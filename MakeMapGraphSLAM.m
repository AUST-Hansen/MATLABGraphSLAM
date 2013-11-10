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
%% Constant parameters
% File to read
filename = 'RandomTrajectories/WallFollowingTraj2.mat';
%filename = 'RandomTrajectories/WideAngleWallFollowingTraj1.mat';
% Number of points to use
% startIdx and endIdx will be used for all data. skip will only be applied to
% range data, to allow a sparser map of observed features.
startIdx = 1; % for now, this must be 1
skip = 1;     % for now, this must be 1 also
rangeSkip = 10;
endIdx = 1000;
stateSize = 6;
% Imagenex matching parameters
ignx_sparsity = 3;
scanmatch_threshold = 30;
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
fprintf('Reson...\n')
[trueBergShape,~] = reconstructTrueCloud(x_obs_t_true(:,startIdx:skip:endIdx),euler_obs_t(:,startIdx:skip:endIdx),x_berg_cm_t(:,startIdx:skip:endIdx),euler_berg_t(:,startIdx:skip:endIdx),sensor,truthrangeData(:,startIdx:skip:endIdx));
%% Clean data of bogus ranges %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if range == -17 or MAX_Range, remove it.
%--------------------------------------------------------------------------
% MeasurementTimeStamps allows subsampling of range data, meaning we can
% integrate high rate pose data but not deal with as many map features if
% we don't want to.
fprintf('Cleaning and formatting data...\n')
[measurementTimestamps, rangeMeasurements, c_i_t] = cleanMeasurements(timeSteps(1:skip:endIdx),sensor,rangeData(:,startIdx:skip:endIdx),useMultibeamRanges,imagenex,imagenexData(:,startIdx:skip:endIdx),useImagenexData,dvl,dvlData(startIdx:skip:endIdx),useDVLRanges,rangeSkip);
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
initialPsiDotGuess = .00;
initialStateEstimate = GraphSLAM_initialize(timeSteps(startIdx:endIdx),inputs,initialPsiDotGuess);
fprintf('Initializing Map Features...\n')
initialMapEstimate = GraphSLAM_initializeMap(initialStateEstimate,rangeMeasurements,c_i_t);
initialFullStateEstimate = [reshape(initialStateEstimate,[],1); reshape(initialMapEstimate,[],1)];
%% Calculate inzformation form of full posterior
fprintf('Linearizing...\n')
[Omega,zeta] = GraphSLAM_linearize(timeSteps(startIdx:endIdx),inputs,measurementTimestamps,imagenex,rangeMeasurements,dvl,dvlData,c_i_t,initialFullStateEstimate);
%% initial testing stuff
%stateHist = reshape(Omega(1:endIdx*6,1:endIdx*6)\zeta(1:endIdx*6),6,[]);
figure(1);
plot(timeSteps,euler_berg_t(3,:));
title('iceberg heading')
state = Omega\zeta;
stateHist = reshape(state(1:6*endIdx),6,[]);
mapEsts = reshape(state(6*endIdx+1:end),3,[]);
figure(2)
plot(-initialStateEstimate(2,:),initialStateEstimate(1,:))
hold on
axis equal
plot(-stateHist(2,:),stateHist(1,:),'r')
plot(xVehBergframe(1,:) - xVehBergframe(1,1)  ,xVehBergframe(2,:) - xVehBergframe(2,1),'g')
scatter(-initialMapEstimate(2,:),initialMapEstimate(1,:),ones(1,size(initialMapEstimate,2)),'b')
scatter(-mapEsts(2,:),mapEsts(1,:),ones(1,size(mapEsts,2)),'k')
legend('initial estimate','after being pushed through information form','truth')
title('estimated path')
%% Reduce graph by marginalizing 
fprintf('Reducing...\n')
[OmegaReduced,zetaReduced] = GraphSLAM_reduce(timeSteps(startIdx:endIdx),stateSize,Omega,zeta);
% state2 = OmegaReduced\zetaReduced;
% stateHist2 = reshape(state2,6,[]);
% figure(1)
% plot(-stateHist2(2,:),stateHist2(1,:),'.k')
% hold on
% axis equal
% title('reduced solution')
%% Initial solve
fprintf('Solving...\n')
[mu, Sigma] = GraphSLAM_solve(OmegaReduced,zetaReduced,Omega,zeta);
% Timeout counter for main loop
itimeout = 0;
MAX_ITER = 2000;
while (itimeout < MAX_ITER) 
    
    % Test correspondences for imagenex scan matching.
    for ii = 1:ignx_sparsity:endIdxtimeSteps
       for jj = ii+ignx_sparsity:ignx_sparsity:endIdxtimeSteps 
           % limit radius of net you cast
           if (norm(mu(1:2,ii) - mu(1:2,jj)) < scanmatch_threshold)
              [validMatch, deltaPose] = correlateImagenex(mu(:,ii),mu(:,jj),imagenexData(ii,:),imagenexData(jj,:),imagenex);
              % if good match, draw link between two
           end
       end
    end
    
    
    
    
   itimeout=itimeout+1 ;
   fprintf('%d iterations completed...\n',itimeout)
end


