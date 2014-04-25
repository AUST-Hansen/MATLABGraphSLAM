clear all; close all; clc;

load wallSegmentPass1.mat
pointCloud1 = pointCloud;
%clear pointCloud
load wallSegmentPass2.mat
pointCloud2 = pointCloud;
clear pointCloud

pointStart = 10000;
pointEnd = 30000;%round(length(pointCloud1));
%pointStart2 = 10000;
%pointEnd2 = 23000; %round(length(pointCloud)/10);
sparsity = 1;


% focus on smaller, subsampled cloud for initial testing
subCloud1 = pointCloud1(:,pointStart:sparsity:pointEnd) + .00*rand(size(pointCloud1(:,pointStart:sparsity:pointEnd)));
%subCloud2 = pointCloud2(:,pointStart2:sparsity:pointEnd2) + .0*rand(size(pointCloud2(:,pointStart2:sparsity:pointEnd2)));

% clean subclouds based on statistical surface characterization
cleanCloud1 = cleanCloud(subCloud1,'kNeighbors',50,'alpha',1);
%cleanCloud2 = cleanCloud(subCloud2,'kNeighbors',500,'alpha',.42);

sprsty = 1;
figure;
scatter3(subCloud1(2,1:sprsty:end),subCloud1(1,1:sprsty:end),-subCloud1(3,1:sprsty:end),3*ones(1,size(subCloud1(1,1:sprsty:end),2)),'b');
axis equal
view(174,40)
figure;
scatter3(cleanCloud1(2,1:sprsty:end),cleanCloud1(1,1:sprsty:end),-cleanCloud1(3,1:sprsty:end),3*ones(1,size(cleanCloud1(1,1:sprsty:end),2)),'r');
axis equal
view(174,40)