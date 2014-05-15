clear all; close all; clc;

if(~exist('M','var'))
	load extractedDataFullWall.mat
end


[timeStamps, IA, IC] = unique(M(:,1));

timeSteps = timeStamps - timeStamps(1);

% first pass: 0 - 5900
% second pass: 6000 - 14000
% third pass: 14400 - 21000

LatLonDepth = M(IA,2:4);
Psi = M(IA,5)*pi/180.;
Speed = M(IA,6)*.3048 ; % I think this is in feet/sec, so convert to m/s

[eastings,northings,zone] = deg2utm(LatLonDepth(:,1),LatLonDepth(:,2));

XYZ = [northings, eastings, LatLonDepth(:,3)];

size(timeSteps)


if(~exist('c_i_t','var'))
    nBeams = 512;
    c_i_t = zeros(nBeams,length(timeSteps));
    rangecounter = 1;
    for ii = 1:length(timeSteps)
        
        range_t = M(IC==ii,7:9)';
        m = length(range_t);
        if (m>nBeams) % duplicate record
            c_i_t(1:m/2,ii) = (rangecounter:rangecounter+(m/2)-1)';
        else
            c_i_t(1:m,ii) = (rangecounter:rangecounter+m-1)';
        end
        rangecounter = rangecounter+m;
        if (mod(ii,100) == 0)
            fprintf('processed %3.1f of %3.1f scans\n',ii,length(timeSteps))
        end
    end
    save(savedmatfile,'M','filename','c_i_t')
else
    nBeams = size(c_i_t,1);
end
%spy(c_i_t);
rangeMeasurements = M(:,7:9)';
%%

figureHandle = 1;
figure(figureHandle);scatter3(XYZ(:,2),XYZ(:,1),-XYZ(:,3));axis equal
title('AUV trajectory in UTM and depth')
axis equal
view(-132,24)
% build point cloud

% read in features
%load extractedFeatures_all
load extractedFeatures5_10meters

% find out most features seen in a timestep
[ifeat, jfeat] = find(correspondences);
maxNumFeatures = max(ifeat);

poseskip = 1;
rangeskip = 5;
startidx = 1;
endidx = length(timeSteps);
pointCloud = zeros(3,length(startidx:poseskip:endidx)*length(1:rangeskip:nBeams));
featureCloud = zeros(size(featurePositions));
pointCloudCounter = 1;
featureCloudCounter = 1;
pointRanges = [];
tic
for ii = startidx:poseskip:endidx % for each pose
    
    % Rotation matrix from vehicle to berg
    position = XYZ(ii,:)';
    heading = Psi(ii);
    i_R_v = Euler2RotMat(0,0,heading);
    
    for jj = 1:rangeskip:nBeams
        measIdx = c_i_t(jj,ii);
        
        if (measIdx > 0)
            pointCloud(:,pointCloudCounter) = position + i_R_v*rangeMeasurements(:,measIdx);
            pointCloudCounter = pointCloudCounter+1;
        end
    end
    
    for jj = 1:maxNumFeatures
        measIdx = correspondences(jj,ii);
        
        if (measIdx > 0)
           featureCloud(:,featureCloudCounter) = position + i_R_v*featurePositions(:,measIdx); 
           featureCloudCounter = featureCloudCounter + 1;
        end
    end
    
end
toc
pointCloud = pointCloud(:,1:pointCloudCounter-1);
%%
pointCloud = pointCloud(:,1:pointCloudCounter-1);
figureHandle = 1;
sprsty = 3;
figure(figureHandle);scatter3(XYZ(1:sprsty:end,2),XYZ(1:sprsty:end,1),-XYZ(1:sprsty:end,3),3*ones(1,size(XYZ(1:sprsty:end,1),2)));axis equal
title('AUV trajectory in UTM and depth')
axis equal
view(-132,24)
hold on;
sparsity = 50;
pointCloudNorm = pointCloud - repmat(mean(pointCloud,2),1,size(pointCloud,2));
featureCloudNorm = featureCloud - repmat(mean(pointCloud,2),1,size(featureCloud,2));

figure(figureHandle)
scatter3(pointCloud(2,1:sparsity:end),pointCloud(1,1:sparsity:end),-pointCloud(3,1:sparsity:end),2*ones(1,size(pointCloud(1,1:sparsity:end))),'g');
scatter3(featureCloud(2,1:sprsty:2248),featureCloud(1,1:sprsty:2244),-featureCloud(3,1:sprsty:2244),'r^');
scatter3(featureCloud(2,2245:sprsty:4828),featureCloud(1,2245:sprsty:4828),-featureCloud(3,2245:sprsty:4828),'k+');
scatter3(featureCloud(2,4829:sprsty:end),featureCloud(1,4829:sprsty:end),-featureCloud(3,4829:sprsty:end),'c*');


