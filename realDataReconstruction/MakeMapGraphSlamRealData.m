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
tEnd = 2000;
processedDataFileName = ['ProcessedData' num2str(tStart) '_' num2str(tEnd) '.mat'];
if (~exist(processedDataFileName,'file')
    if (~exist(savedfilename,'file'))
        ExtractedData = extractData(filename);
        save(savedfilename,'ExtractedData');
    else
        load(savedfilename)
    end
    %% Process Data
    %
    %   Inputs:
    %       ExtractedData
    %       tStart
    %       tEnd
    %
    ProcessedData = processData(ExtractedData,tStart,tEnd);
    save(processedDataFileName,ProcessedData)
    clear ExtractedData
else
    load(processedDataFileName)
end

% Extract position

keyboard

% OtherData contains data like rangeMeasurements, which are used often, by
% many functions, but not changed during the course of GraphSLAM.
%% Initialize SLAM
[SLAMdata, OtherData] = GraphSLAMInitializeRD(ProcessedData); % RD = real data
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

