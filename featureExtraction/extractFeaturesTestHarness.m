clear all; close all; clc;

load wallSegmentPass1.mat
pointCloud1 = pointCloud;
%clear pointCloud
load wallSegmentPass2.mat
pointCloud2 = pointCloud;
clear pointCloud
DEBUG = true;

pointStart = 10000;
pointEnd = 24000; %round(length(pointCloud)/10);
pointStart2 = 8000;
pointEnd2 = 25000; %round(length(pointCloud)/10);
% pointStart = 17000;
% pointEnd = 22000; %round(length(pointCloud)/10);
% pointStart2 = 19000;
% pointEnd2 = 23000; %round(length(pointCloud)/10);
sparsity = 6;


% focus on smaller, subsampled cloud for initial testing
subCloud1 = pointCloud1(:,pointStart:sparsity:pointEnd) + .00*rand(size(pointCloud1(:,pointStart:sparsity:pointEnd)));
subCloud2 = pointCloud2(:,pointStart2:sparsity:pointEnd2) + .0*rand(size(pointCloud2(:,pointStart2:sparsity:pointEnd2)));
subCloud1 = subCloud1 - repmat(mean(subCloud1,2),1,size(subCloud1,2));
subCloud2 = subCloud2 - repmat(mean(subCloud2,2),1,size(subCloud2,2));
sprsty = 1;

% clean subclouds based on statistical surface characterization
subCloud1 = cleanCloud(subCloud1,'kNeighbors',40,'alpha',1);
subCloud2 = cleanCloud(subCloud2,'kNeighbors',40,'alpha',1);

if (DEBUG)
    figure
    offset = 0*ones(size(subCloud2));
    offset(3,:) = 0;
    showCloud2 = subCloud2+offset;
    
    axis equal
    hold on
    scatter3(subCloud1(2,1:sprsty:end),subCloud1(1,1:sprsty:end),-subCloud1(3,1:sprsty:end),3*ones(1,size(subCloud1(1,1:sprsty:end),2)),'c');
    scatter3(showCloud2(2,1:sprsty:end),showCloud2(1,1:sprsty:end),-showCloud2(3,1:sprsty:end),3*ones(1,size(subCloud2(1,1:sprsty:end),2)),'g^');
    view(-165,6)
end
%subCloud1 = subCloud1 - repmat(mean(subCloud1,2),1,size(subCloud1,2));
%subCloud2 = subCloud2 - repmat(mean(subCloud2,2),1,size(subCloud2,2));
%%
tic
[measurements1, descriptors1] = extractFeaturesFromSubmap(subCloud1,'Verbose',false,'minRadius',3,'maxRadius',7,'numRadii',4,'Sparsity',5,'PFHbins',4);
toc
[measurements2, descriptors2] = extractFeaturesFromSubmap(subCloud2,'Verbose',false,'minRadius',3,'maxRadius',7,'numRadii',4,'Sparsity',5,'PFHbins',4);

%%
map = figure();
sprsty = 1;

%scatter3(measurements1(2,:),measurements1(1,:),-measurements1(3,:),'k.');
for ii = 1:length(measurements1)
    %text(measurements1(2,ii),measurements1(1,ii),-measurements1 (3,ii),num2str(ii))
end
offset = 40*ones(size(subCloud2));
offset(3,:) = -20;

showCloud2 = subCloud2+offset;
showMeas2 = measurements2+offset(:,1:length(measurements2));

axis equal
hold on
scatter3(subCloud1(2,1:sprsty:end),subCloud1(1,1:sprsty:end),-subCloud1(3,1:sprsty:end),3*ones(1,size(subCloud1(1,1:sprsty:end),2)),'c');
scatter3(measurements1(2,:),measurements1(1,:),-measurements1(3,:),'ko');
scatter3(showCloud2(2,1:sprsty:end),showCloud2(1,1:sprsty:end),-showCloud2(3,1:sprsty:end),3*ones(1,size(subCloud2(1,1:sprsty:end),2)),'g^');
scatter3(showMeas2(2,:),showMeas2(1,:),-showMeas2(3,:),'rx');

for ii = 1:length(showMeas2)
    %text(showMeas2(2,ii),showMeas2(1,ii),-showMeas2(3,ii),num2str(ii))
end
axis equal
view(-165,6)
%%
matches = zeros(1,size(descriptors2,2));
absthresh = .9;
relthresh = 1.;
maxdist = max( norm(max(measurements1,[],2)-min(measurements1,[],2)),norm(max(measurements2,[],2)-min(measurements2,[],2)));
for ii = 1:size(descriptors1,2)
    bestidx = 0;
    best = -1;
    nextbest = -1;
    for jj = 1:size(descriptors2,2)
        % inner product
        ip = descriptors1(:,ii)'*descriptors2(:,jj); 
        if(false)
            plot(descriptors1(:,ii),'r');
            hold on;
            plot(descriptors2(:,jj),'b');
            drawnow;
            hold off;
            pause(.01);
        end
        if ip>best
            nextbest = best;
            best = ip;
            bestidx = jj;
        end
    end
    [best nextbest/best]
    if and(best > absthresh, (best)/(nextbest+eps) > relthresh)
        % declare match
        matches(ii) = bestidx;
    end
end
figure(map);
for ii = 1:min(size(descriptors1,2),size(descriptors2,2))
    if matches(ii) > 0 % goodmatch
        hold on;
        lines = [measurements1(:,ii)' ; showMeas2(:,matches(ii))'];
        [ii matches(ii) descriptors1(:,ii)'*descriptors2(:,matches(ii))]
        plot3(lines(:,2),lines(:,1),-lines(:,3),'r')
        hold on;
        %drawnow()
        %pause(.2)
    end
end
matches;
% 
% %% RANSAC
% [R,T,inliers] = ransacPointCloud(measurements1,showMeas2,matches);
% 
% q4 = R*showMeas2 + repmat(T,1,length(showMeas2));
% figure(map);
% 
% scatter3(q4(2,:),q4(1,:),-q4(3,:),'k^');
%%
%SC1 = subCloud1 - repmat(mean(subCloud1,2),1,size(subCloud1,2));
[TR, TT] = icp(subCloud1,subCloud2,'WorstRejection',.5);

plotsparsity = 1;
%subCloud1 = pointCloud1(:,pointStart:plotsparsity:pointEnd);
%subCloud2 = pointCloud2(:,pointStart:plotsparsity:pointEnd);
figure;
axis equal;
hold on
scatter3(subCloud1(2,1:sprsty:end),subCloud1(1,1:sprsty:end),-subCloud1(3,1:sprsty:end),3*ones(1,size(subCloud1(1,1:sprsty:end),2)),'c');
scatter3(measurements1(2,:),measurements1(1,:),-measurements1(3,:),'ko');
axis equal;
hold on
dispCloud2 = TR*subCloud2 + repmat(TT,1,size(subCloud2,2));
measCloud2 = TR*measurements2 + repmat(TT,1,size(measurements2,2));
scatter3(dispCloud2(2,1:sprsty:end),dispCloud2(1,1:sprsty:end),-dispCloud2(3,1:sprsty:end),3*ones(1,size(dispCloud2(1,1:sprsty:end),2)),'g^');
scatter3(measCloud2(2,:),measCloud2(1,:),-measCloud2(3,:),'r+');
view(178,20)
