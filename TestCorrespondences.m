%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    TestCorrespondences.m
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
poseSkip = 1;
endIdx = 1000;
stateSize = 6;
% Imagenex matching parameters
ignx_sparsity = 3;
scanmatch_threshold = 100; % only look for matches within this many meters of pose difference
scanmatch_RMStolerance = 2;
scanmatch_probThreshold = .5; 
% Use imagenex, multibeam and DVL ranges in observations?
useDVLRanges = false;
useMultibeamRanges = false;
useImagenexData = true;

%% Read in Data
fprintf('Reading in Data...\n')
addpath ..
load(filename)
%% Clean data of bogus ranges %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if range == -17 or MAX_Range, remove it.
%--------------------------------------------------------------------------
% MeasurementTimeStamps allows subsampling of range data, meaning we can
% integrate high rate pose data but not deal with as many map features if
% we don't want to.
fprintf('Cleaning and formatting data...\n')
[measurementTimestamps, rangeMeasurements, c_i_t] = cleanMeasurements(timeSteps(1:skip:endIdx),sensor,rangeData(:,startIdx:skip:endIdx),useMultibeamRanges,imagenex,imagenexData(:,startIdx:skip:endIdx),useImagenexData,dvl,dvlData(startIdx:skip:endIdx),useDVLRanges,rangeSkip,poseSkip);
%--------------------------------------------------------------------------
%% Set up icp stuff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
residSurface = zeros(length(measurementTimestamps));
for jj = 1:size(measurementTimestamps,2)-1
pose1FeatureIndex = find (c_i_t(:,jj)~=-17);
featuresSeenAtPose1 = c_i_t(pose1FeatureIndex,jj);
pointCloud1 = rangeMeasurements(:,featuresSeenAtPose1);
resids = 0*measurementTimestamps;

for ii = jj+1:size(measurementTimestamps,2)
    % build poses
    pose2FeatureIndex = find (c_i_t(:,ii)~=-17);
    featuresSeenAtPose2 = c_i_t(pose2FeatureIndex,ii);
    pointCloud2 = rangeMeasurements(:,featuresSeenAtPose2);
    % do icp
    [R,T,ERR] = icp(pointCloud1,pointCloud2);
    
    resids(ii) = ERR(end);
    if(mod(ii,50) == 0)
       %fprintf('%d of %d poses\n',ii,size(measurementTimestamps,2)) ;
    end
end
fprintf('%d\n',jj)
residSurface(jj,:) = resids;
end
%%
plot(residSurface(1,:))
figure;
surf(residSurface)
lcResids = resids;
lcResids(1:50) = 1000;
%% look for loop closure
[y,i]=min(lcResids);
pose2FeatureIndex = find (c_i_t(:,i)~=-17);
featuresSeenAtPose2 = c_i_t(pose2FeatureIndex,i);
pointCloud2 = rangeMeasurements(:,featuresSeenAtPose2);
% do icp
[R,T,ERR] = icp(pointCloud1,pointCloud2);
pNew = R*pointCloud2 + repmat(T,1,size(pointCloud2,2));
figure(9)
hold on; axis equal;
for qq = 1:length(featuresSeenAtPose1)
    scatter(pointCloud1(1,qq),pointCloud1(2,qq),'r')
    text(pointCloud1(1,qq),pointCloud1(2,qq)+1,num2str(featuresSeenAtPose1(qq)));
end
for qq= 1:length(featuresSeenAtPose2)
    scatter(pNew(1,qq),pNew(2,qq),'b')
    text(pNew(1,qq),pNew(2,qq)-1,num2str(featuresSeenAtPose2(qq)));
end