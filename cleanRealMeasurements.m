function [measurementTimestamps, rangeMeasurements, correspondences,meas_ind,multibeamMeasurements,correspondences_mb] = cleanRealMeasurements(filename)
%clear all; clc;

savedmatfile = 'extractedDataFullWall.mat';

if(~exist(savedmatfile,'file'))
    M = csvread(filename,1,0);
    save(savedmatfile,'M','filename')
else
    load(savedmatfile)
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

figure;scatter3(XYZ(:,1),XYZ(:,2),XYZ(:,3));axis equal
title('AUV trajectory in UTM and depth')

% Plot data

startindices = [3000 10000];
endindices = [4000 11000];
for qq = 1:2
    poseskip = 2;
    rangeskip = 5;
    startidx = startindices(qq); %3000;
    endidx = endindices(qq); %4000;%length(timeSteps);
    pointCloud = zeros(3,length(startidx:poseskip:endidx)*length(1:rangeskip:nBeams));

    pointCloudCounter = 1;
    pointRanges = [];
    for ii = startidx:poseskip:endidx % for each pose
        
        % Rotation matrix from vehicle to berg
        position = XYZ(ii,:)';
        heading = Psi(ii);
        i_R_v = Euler2RotMat(0,0,heading);
        
        for jj = 1:rangeskip:nBeams
            measIdx = c_i_t(jj,ii);
            
            if (measIdx > 0);
                pointCloud(:,pointCloudCounter) = position + i_R_v*rangeMeasurements(:,measIdx);
                pointCloudCounter = pointCloudCounter+1;
            end
        end
        
        if ((ii<endidx - 2*poseskip) && timeSteps(ii+poseskip) - timeSteps(ii) > 400)
            pointRanges = [pointRanges pointCloudCounter]
        end
        
    end
    
    % window different passes
    discontinuityIndices = find(diff(timeSteps) > 100);
    nPasses = length(pointRanges) + 1;
    keyboard;
    pointCloud = pointCloud(:,1:pointCloudCounter-1);
    save(['wallSegmentPass' num2str(qq) '.mat'],'pointCloud');    
    tNavStart = 1;
    tNavEnd = length(timeSteps);%discontinuityIndices(1);
    pointStart = 1;
    if (isempty(pointRanges))
        pointEnd = pointCloudCounter-1;
    else
        pointEnd = pointRanges(1);
    end
    colorz = 'rbk';
    colorzz = 'mcg';
    for iPass = 1:nPasses
        scatter3(pointCloud(2,pointStart:pointEnd),pointCloud(1,pointStart:pointEnd),-pointCloud(3,pointStart:pointEnd),2*ones(1,size(pointCloud(1,pointStart:pointEnd))),colorz(qq));
        hold on;
        
        pointStart = pointEnd+1;
        
        if (iPass < length(pointRanges))
            pointEnd = pointRanges(iPass + 1);
            
        else
            pointEnd = length(pointCloud);
            
        end
        iPass
    end
end    
    for iPass = 1:nPasses
        
        scatter3(XYZ(tNavStart:tNavEnd,2),XYZ(tNavStart:tNavEnd,1),-XYZ(tNavStart:tNavEnd,3),colorzz(iPass))
        
        tNavStart = tNavEnd+1;
        if (iPass < length(discontinuityIndices))
            
            tNavEnd = min(discontinuityIndices(iPass + 1),endidx);
        else
            
            tNavEnd = endIdx;
        end
        iPass
    end

view(122,15)
axis equal




%    scatter3(-pointCloud(1,1:floor(end/2)),-pointCloud(2,1:floor(end/2)),-pointCloud(3,1:floor(end/2)),2*ones(1,size(pointCloud(1,1:floor(end/2)))));
%     hold on;
%    scatter3(-pointCloud(1,floor(end/2)+1:end),-pointCloud(2,floor(end/2)+1:end),-pointCloud(3,floor(end/2)+1:end),2*ones(1,size(pointCloud(1,floor(end/2)+1:end))),'r');
%scatter3(-XYZ(1:floor(end/2),1),-XYZ(1:floor(end/2),2),-XYZ(1:floor(end/2),3),'c')
%scatter3(-XYZ(floor(end/2)+1:end,1),-XYZ(floor(end/2)+1:end,2),-XYZ(floor(end/2)+1:end,3),'k')

