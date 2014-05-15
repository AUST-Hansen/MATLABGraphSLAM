function [featureMeasurements, featureDescriptors, featureCorrespondences] = ExtractPseudomeasurements(M,c_i_t,spacing,startTime,endTime,varargin)

%% Use the inputParser class to validate input arguments.
%inp = inputParser;

%inp.addRequired('M');

%inp.addRequired('c_i_t');

%inp.addRequired('spacing', @(x)isscalar(x) && x > 0);

%inp.addRequired('startTime', @(x)isscalar(x) && x >= 0);

%inp.addRequired('endTime', @(x)isscalar(x) && x > 0);

%inp.addOptional('poseskip', 1, @(x)isscalar(x) && x > 0); % multiplier on radius

%inp.addOptional('rangeskip', 1, @(x)isscalar(x) && x > 0);

%validMinimize = {'point','plane','lmapoint'};
%inp.addParamValue('Minimize', 'point', @(x)any(strcmpi(x,validMinimize)));
%inp.addParamValue('Triangulation', [], @(x)isreal(x) && size(x,2) == 3);

%inp.addParamValue('Verbose', false, @(x)islogical(x));

%inp.parse(varargin{:});
%arg = inp.Results;
%clear('inp');
arg.startTime = startTime;
arg.spacing = spacing;
arg.endTime = endTime;
arg.rangeskip = 2;
arg.poseskip = 1;

% if(~exist('extractedDataFullWall.mat','file'))
%     M = csvread(filename,1,0);
%     save('extractedDataFullWall.mat','M','filename')
% else
%     load extractedDataFullWall.mat
% end

[timeStamps, IA, IC] = unique(M(:,1));

timeSteps = timeStamps - timeStamps(1);

if startTime < timeSteps(1)
    fprintf('Start time out of range. First time on record = %d\n',timeSteps(1));
    fprintf('setting Start time to %d...\n',timeSteps(1));
    startTime = timeSteps(1);
end
if endTime > timeSteps(end)
    fprintf('End time out of range. Last time on record = %d\n',timeSteps(1));
    fprintf('setting End time to %d...\n',timeSteps(1));
    endTime = timeSteps(end);    
end

startIndex = find(timeSteps>=startTime,1);
if isempty(startIndex)
    fprintf('Warning 1\n');
    startIndex = 1;
end
endIndex = find(timeSteps<=endTime,1,'last');
if isempty(endIndex)
    fprintf('Warning 2\n');
    endIndex = length(timeSteps);
end

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
    save('extractedDataWall.mat','M','filename','c_i_t')
else
    nBeams = size(c_i_t,1);
end
%spy(c_i_t);
%rangeMeasurements = M(:,7:9)';
rangeMeasurements = M( (M(:,1) - M(1,1) >=arg.startTime & M(:,1) - M(1,1) <=arg.endTime) ,7:9)';

subcit = c_i_t(:,(timeSteps >=arg.startTime & timeSteps <=arg.endTime));
firstmeas = min(subcit(subcit~=0));

%%
%figure;scatter3(XYZ(:,2),XYZ(:,1),-XYZ(:,3));axis equal
%title('AUV trajectory in UTM and depth')

% c_i_t for features

featureCorrespondences = spalloc(nBeams,length(timeSteps),round((length(timeSteps)/spacing)*200));
featureMeasurements = [];
featureDescriptors = [];
featureCounter = 0;

for iRefFrame = startIndex+spacing:spacing:endIndex-spacing
    %iRefFrame
    timeSteps(iRefFrame)
    pointCloud = zeros(3,2*spacing+10*length(1:arg.rangeskip:nBeams));
    pointCloudCounter = 1;
    %pointRanges = [];
    for ii = iRefFrame-spacing:arg.poseskip:iRefFrame+spacing
        % Rotation matrix from vehicle to berg
        position = XYZ(ii,:)';
        heading = Psi(ii);
        % IT WOULD BE REALLY NICE TO HAVE PITCH AND ROLL!!!!!!!!!!!
        i_R_v = Euler2RotMat(0,0,heading);
        
        for jj = 1:arg.rangeskip:nBeams
            measIdx = c_i_t(jj,ii)-firstmeas+1;
            
            if (measIdx > 0);
                pointCloud(:,pointCloudCounter) = position + i_R_v*rangeMeasurements(:,measIdx);
                pointCloudCounter = pointCloudCounter+1;
            end
        end
        
    
    
    end % ii
    % Trim point cloud
    pointCloudWorldFrame = pointCloud(:,1:pointCloudCounter-1);
    % put into vehicle frame
    vRi = Euler2RotMat(0,0,Psi(iRefFrame))';
    pointCloudVehicleFrame = vRi*(pointCloudWorldFrame - repmat(XYZ(iRefFrame,:)',1,size(pointCloudWorldFrame,2)));
    
    
    % Extract pseudomeasurements
    [measurements_t, descriptors_t] = extractFeaturesFromSubmap(pointCloudVehicleFrame(:,1:20:end),'Verbose',false,'minRadius',5,'maxRadius',10,'numRadii',4,'Sparsity',5,'PFHbins',4);
    featureMeasurements = [featureMeasurements measurements_t];
    featureDescriptors = [featureDescriptors descriptors_t];
    featureCorrespondences(1:size(measurements_t,2),iRefFrame) = (featureCounter+1:featureCounter+size(measurements_t,2))';
    featureCounter = featureCounter + size(measurements_t,2);

    if(false)
        figure(3);scatter3(XYZ(:,2),XYZ(:,1),-XYZ(:,3));axis equal
        hold on;
        title('AUV trajectory in UTM and depth')
        scatter3(pointCloudWorldFrame(2,:),pointCloudWorldFrame(1,:),-pointCloudWorldFrame(3,:),ones(1,size(pointCloudWorldFrame,2)),'r')
        axis equal
        view(-106,26)
        hold off
        figure(4)
        scatter3(pointCloudVehicleFrame(2,:),pointCloudVehicleFrame(1,:),-pointCloudVehicleFrame(3,:),ones(1,size(pointCloudWorldFrame,2)),'r')
        axis equal;
        view(-106,26)
        drawnow()
        pause(.1)
        hold off;
    end
end % iRefFrame
