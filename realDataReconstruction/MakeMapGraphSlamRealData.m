clear all; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   MakeMapGraphSlamRealData.m
%
%   Takes in real data, performs GraphSLAM 
%   State = [x y psi u v psidotbias]';
%   Costate = [z phi theta]' -- Assumed perfect knowledge
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Constant parameters
% File to read
%filename = 'RandomTrajectories/WallFollowingTraj2.mat';
trajname = 'fullWallextractedData';
location = '.';
savedfilename = [location '/' trajname '.mat'];
filename = [location '/' trajname '.csv'];
% GraphSLAM control parameters
GSparams.MAX_ITER = 10;
GSparams.iteration = 1;

%
%% Extract data
tStart = 0;
tEnd = 5200;
submaphalfwidth = 400;
processedDataFileName = ['ProcessedData' num2str(tStart) '_' num2str(tEnd) '.mat'];
featureDataFileName = ['PDfeatures' num2str(tStart) '_' num2str(tEnd) '.mat'];
if (exist(featureDataFileName))
    load featureDataFileName
else
    if (~exist(processedDataFileName,'file'))
        fprintf('No processed data file found: Looking for raw data file...\n');
        if (~exist(savedfilename,'file'))
            fprintf('Extracting raw data from csv...\n');
            ExtractedData = extractData(filename);
            fprintf('saving to mat file\n');
            save(savedfilename,'ExtractedData');
        else
            fprintf('MAT file found: Loading raw data from file...\n');
            load(savedfilename)
        end
        %% Process Data
        %
        %   Inputs:
        %       ExtractedData
        %       tStart
        %       tEnd
        %
        fprintf('Processing raw data...\n');
        ProcessedData = processData(ExtractedData,tStart,tEnd,submaphalfwidth);
        fprintf('Saving processed data...\n');
        save(processedDataFileName,'ProcessedData')
        clear ExtractedData
    else
        fprintf('Processed data found! Loading from file...\n');
        load(processedDataFileName)
    end
    PDfeatures = sweepSonarFeatures(ProcessedData);
    save(featureDataFileName,'PDfeatures')
    clear ProcessedData
end

% Extract position



% OtherData contains data like rangeMeasurements, which are used often, by
% many functions, but not changed during the course of GraphSLAM.
%% Initialize SLAM
[SLAMdata, OtherData] = GraphSLAMInitializeRD(PDfeatures); % RD = real data

SLAMdata = GraphSLAMLinearizeRD(SLAMdata,OtherData);
SLAMdata = GraphSLAMReduceRD(SLAMdata,OtherData);
SLAMdata = GRAPHSLAMSolveRD(SLAMdata,OtherData);
%% run SLAM
while (GSparams.iteration < GSparams.MAX_ITER)
    % Track innovations
    LastSLAMdata = SLAMdata;
    % Test for new correspondences
    SLAMdata = GraphSLAMCorrespondenceTest(SLAMdata,OtherData);
    % Run slam
    SLAMdata = GraphSLAMLinearizeRD(SLAMdata,OtherData);
    SLAMdata = GraphSLAMReduceRD(SLAMdata,OtherData);
    SLAMdata = GRAPHSLAMSolveRD(SLAMdata,OtherData);
    
    
   GSparams.iteration = GSparams.iteration + 1; 
end

