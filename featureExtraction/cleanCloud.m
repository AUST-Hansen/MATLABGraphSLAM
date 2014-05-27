function [cleanedCloud] = cleanCloud(cloud,varargin)

    % cloud must be expressed in the frame of reference of the vehicle at
    % the reference time, with the origin of the cloud at the vehicle
    % location.
    
%% Use the inputParser class to validate input arguments.
inp = inputParser;

inp.addRequired('cloud', @(x)isreal(x) && size(x,1) == 3);

inp.addOptional('kNeighbors', 30.0, @(x)isscalar(x) && x > 0);

inp.addOptional('alpha', 1.0, @(x)isscalar(x) && x > 0); % multiplier on radius

inp.addOptional('Method','Normal');

%validMinimize = {'point','plane','lmapoint'};
%inp.addParamValue('Minimize', 'point', @(x)any(strcmpi(x,validMinimize)));
%inp.addParamValue('Triangulation', [], @(x)isreal(x) && size(x,2) == 3);

inp.addParamValue('Verbose', false, @(x)islogical(x));


inp.parse(cloud,varargin{:});
arg = inp.Results;
clear('inp');


cleanedCloud = zeros(size(cloud));
cleanCounter = 0;

if strcmp(arg.Method,'RegionGrow')
    % clean cloud by growing regions
    numRegions = 0;
    maxRegions = 10;
    sizeRegions = zeros(1,maxRegions);
    
else
    
    % for each point in input cloud find k nearest neighbors
    [idx dK] = knnsearch(cloud',cloud','K',arg.kNeighbors,'Distance','chebychev');
    dists = mean(dK,2)
    
    cov = std(dists)
    
    
    
    
    inliers = dists > mean(dists) - arg.alpha*cov & dists < mean(dists) + arg.alpha*cov;
    
    cleanedCloud = cloud(:,inliers);
    
end
if arg.Verbose
    keyboard
end