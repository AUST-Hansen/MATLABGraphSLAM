function [measurements, descriptors] = extractFeaturesFromSubmap(cloud,varargin)

    % cloud must be expressed in the frame of reference of the vehicle at
    % the reference time, with the origin of the cloud at the vehicle
    % location.
    
%% Use the inputParser class to validate input arguments.
inp = inputParser;

inp.addRequired('cloud', @(x)isreal(x) && size(x,1) == 3);

inp.addOptional('Radius', 5.0, @(x)isscalar(x) && x > 0);

inp.addOptional('MinNeighbors', 15, @(x)isscalar(x) && x > 0); % multiplier on radius

inp.addOptional('Sparsity', 1, @(x)isinteger(x) && x > 0);

inp.addOptional('MinCurve', .06, @(x)isscalar(x) && x > 0);

inp.addOptional('MaxCurve', .2, @(x)isscalar(x) && x > 0);

inp.addOptional('PFHbins',5,@(x)isinteger && x > 0);

%validMinimize = {'point','plane','lmapoint'};
%inp.addParamValue('Minimize', 'point', @(x)any(strcmpi(x,validMinimize)));
%inp.addParamValue('Triangulation', [], @(x)isreal(x) && size(x,2) == 3);

inp.addParamValue('Verbose', false, @(x)islogical(x));


inp.parse(cloud,varargin{:});
arg = inp.Results;
clear('inp');

%%
normals = zeros(size(cloud));
curvatures = zeros(2,length(cloud));

% for each point in
neighborIndices = rangesearch(cloud',cloud',arg.Radius);
for ii = 1:length(cloud)
    
    % weed out borders and spurious points
    if(length(neighborIndices{ii})<arg.MinNeighbors)
        continue;
    end
    
    % make neighbor matrix
    P = cloud(:,neighborIndices{ii});
    Po = mean(P,2);
    Preg = P - repmat(Po,1,size(P,2));
    M = Preg*Preg';
    [u s v] = svd(M);
    normal = v(:,3);
    % Make sure normal is pointed back at vehicle. This check will only get weird
    % with really funky geometries, and should never really cause a
    % problem, based on how the cloud is defined (in vehicle ref frame).
    if (normal'*Po > 0)
        normal = -normal;
    end
    normals(:,ii) = normal;
    curvatures(1,ii) = s(3,3)/(s(3,3)+s(2,2));
    curvatures(2,ii) = s(3,3)/(s(3,3)+s(1,1));
    
    if(arg.Verbose)
        scatter3(Preg(2,:),Preg(1,:),-Preg(3,:),'r');
        hold on;
        axis equal
        quiver3(0,0,0,normal(2),normal(1),-normal(3));
        keyboard
    end
    
end

% have to do the above so that I have all the normals for calculating Point
% Feature Histograms

for ii = 1:length(cloud)
    if and(and(curvatures(1,ii) > arg.MinCurve, true), and(curvatures(2,ii) > arg.MinCurve, true ))   
        % good keypoint. Now describe it
        P = cloud(:,neighborIndices{ii});
        patchNormals = normals(:,neighborIndices{ii});
        descriptor = zeros(arg.PFHbins^3,1);
        for jj = 1:length(P)-1
            for kk = jj+1:length(P)
                n1 = patchNormals(:,jj);
                n2 = patchNormals(:,kk);
                p1 = P(:,jj);
                p2 = P(:,kk);
                % start calculating point feature histogram
                unitDP = (p2-p1)/norm(p2-p1);
                u = n1;
                v = cross(u,unitDP);
                w = cross(u,v);
                alpha = v'*n2; % -1:1
                phi = u'*unitDP; % -1:1
                theta = atan2(w'*n2,u'*n2); %-pi:pi
                % now add this to bins
            end
        end
    end
    
end
