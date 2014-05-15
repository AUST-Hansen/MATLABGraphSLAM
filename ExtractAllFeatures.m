clear all; close all; clc;

addpath featureExtraction

tMax = 500;%8982;
range = 200;
timeBtwRefs = 30; % in num of timesteps. Essentiall divide by 3 for dtime
% this script breaks full wall data into chunks and processes them to
% extract features

featuresSoFar = 0;
featurePositions = [];
featureDescriptors = [];

load extractedDataFullWall.mat

for t = 100:range:tMax
    
    [meas,desc,corrs] = ExtractPseudomeasurements(M,c_i_t,timeBtwRefs,t,t+range);
    
    % add to bin
    featurePositions = [featurePositions meas];
    featureDescriptors = [featureDescriptors desc];
    
    if exist('correspondences','var')
        corrs(corrs>0) = corrs(corrs>0) + featuresSoFar; 
        correspondences = correspondences + corrs;
    else
        correspondences = corrs;
    end
    
    featuresSoFar = featuresSoFar + size(meas,2);
    
end
%%
save('extractedFeatures5_10meters.mat','featurePositions','featureDescriptors','correspondences')