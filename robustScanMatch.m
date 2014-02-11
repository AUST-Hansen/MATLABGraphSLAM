function [TR, TT, ER] = robustScanMatch(pointCloud1,pointCloud2,varargin)

% Default parameter values
kMax = 10;
minPoints = 3;
minAddInliers = 2;
inlierRadius = 1;
planeFitParam = 0.3; % determined empirically

smoothingFactor = 7;
worstRejection = 0.2;
binCount = 5;

% Parse varargin
nVarargs = length(varargin);
if(~isempty(varargin))
    for k = 1:nVarargs
        switch(varargin{k})
            case{'SmoothingFactor'}
                smoothingFactor = varargin{k+1};
            case{'MinPoints'}
                minPoints = varargin{k+1};
            case{'WorstRejection'}
                worstRejection = varargin{k+1};
            case{'InlierRadius'}
                inlierRadius = varargin{k+1};
            case{'PlaneFitParam'}
                planeFitParam = varargin{k+1};
        end %switch
    end %for
end %if

% get length of data
lenP = length(pointCloud1);
lenQ = length(pointCloud2);
% smooth data
smoothingFactor = 10;
smooth1 = [smooth(pointCloud1(1,:),smoothingFactor)';smooth(pointCloud1(2,:),smoothingFactor)';smooth(pointCloud1(3,:),smoothingFactor)'];
smooth2 = [smooth(pointCloud2(1,:),smoothingFactor)';smooth(pointCloud2(2,:),smoothingFactor)';smooth(pointCloud2(3,:),smoothingFactor)'];
%%%%%%%%%%%%%%%%%%%%
% plot for DEBUG
%%%%%%%%%%%%%%%%%%%%
figure;
scatter(pointCloud1(1,:),pointCloud1(2,:),'.b');
hold on;
plot(smooth1(1,:),smooth1(2,:),'-b')
scatter(pointCloud2(1,:),pointCloud2(2,:),'.r');
%hold on;
plot(smooth2(1,:),smooth2(2,:),'-r')
%%%%%%%%%%%%%%%%%%%%
% calculate normals for histogramming
smooth1Prime =  diff(smooth1,1,2);
smooth2Prime = diff(smooth2,1,2);
normalAngles1 = atan2(-smooth1Prime(1,:),smooth1Prime(2,:));
normalAngles2 = atan2(-smooth2Prime(1,:),smooth2Prime(2,:));
normalAngles1 = [normalAngles1 normalAngles1(end)];
normalAngles2 = [normalAngles2 normalAngles2(end)];
%%%%%%%%%%%%%%%%%%%%
% plot for DEBUG
%%%%%%%%%%%%%%%%%%%%
quiver(smooth1(1,1:end-1),smooth1(2,1:end-1),-smooth1Prime(2,:),smooth1Prime(1,:))
axis equal
figure
plot(normalAngles1)
hold on
plot(normalAngles2,'r')
%%%%%%%%%%%%%%%%%%%%
% bin into histogram
range1 = linspace(min(normalAngles1),max(normalAngles1),binCount);
range2 = linspace(min(normalAngles2),max(normalAngles2),binCount);
[bc1, id1] = histc(normalAngles1,range1);
[bc2, id2] = histc(normalAngles2,range2);
% reject outliers
nGoodBins1 = sum(bc1>2);
nGoodBins2 = sum(bc2>2);
minPcount1 = min(bc1(bc1>2));
minPcount2 = min(bc2(bc2>2));
minPcount = min(minPcount1,minPcount2);
pc1 = [];
pc2 = [];
for ii = 1:binCount
    if (bc1(ii) > 2)
        %points = pointCloud1(:,id1 == ii);
        points = smooth1(:,id1 == ii);
        y = randsample(bc1(ii),minPcount);
        pc1 = [pc1 points(:,y)];
    end
        if (bc2(ii) > 2)
        %points = pointCloud2(:,id2 == ii);
        points = smooth2(:,id2 == ii);
        y = randsample(bc2(ii),minPcount);
        pc2 = [pc2 points(:,y)];
    end
end
figure;
scatter(pc1(1,:),pc1(2,:),'b');
hold on
scatter(pc2(1,:),pc2(2,:),'r');
% Sample from bins

% Do ICP
%[TR,TT,ER] = icp(pointCloud1,pointCloud2,20,'WorstRejection',worstRejection,'twoDee',twoDee,'Weight',@weighting);
%[TR,TT,ER] = icp(smooth1,smooth2,20,'WorstRejection',worstRejection);
[TR,TT,ER] = icp(pc1,pc2,10,'WorstRejection',worstRejection);
%cloud2 = TR*pointCloud2 + repmat(TT,1,lenQ);
cloud2 = TR*smooth2 + repmat(TT,1,lenQ);

figure;
offset= 5
%scatter(pointCloud1(1,1:end),pointCloud1(2,1:end),'.b');
scatter(smooth1(1,1:end),smooth1(2,1:end),'.b');
hold on;
scatter(cloud2(1,:),cloud2(2,:),'.r');

[idxKNN, dKNN]=knnsearch(smooth1',cloud2')
%[idxKNN, dKNN]=knnsearch(pointCloud1',cloud2')
[idxSubSample, dSubSample] = knnsearch(pc1',transpose(TR*pc2 + repmat(TT,1,size(pc2,2))))

addInliers = sum(dKNN<inlierRadius) - sum(dSubSample < inlierRadius)
dKNN = min(dKNN,inlierRadius);
%figure; plot(dKNN)
norm(dKNN)

ER;




    function weights = weighting(indices)
        indices
        size(w1)
        size(w2)
        w1
        w2
        weights = w1(indices).*w2;
        

