function [TR, TT, ER] = ransacICP2d(pointCloud1,pointCloud2,varargin)


TR = [];
TT = [];
ER = [];

sz1 = size(pointCloud1,2);
sz2 = size(pointCloud2,2);
cloudsize = min(sz1,sz2);

% Default parameter values
kMax = 10;
minPoints = 3;
minAddInliers = 2;
inlierRadius = 1;
planeFitParam = 0.3; % determined empirically

% Parse varargin
nVarargs = length(varargin);
if(~isempty(varargin))
    for k = 1:nVarargs
        switch(varargin{k})
            case{'MaxIter'}
                kMax = varargin{k+1};
            case{'MinPoints'}
                minPoints = varargin{k+1};
            case{'MinAddInliers'}
                minAddInliers = varargin{k+1};
            case{'InlierRadius'}
                inlierRadius = varargin{k+1};
            case{'PlaneFitParam'}
                planeFitParam = varargin{k+1};
        end %switch
    end %for
end %if

kMax
minAddInliers
inlierRadius
minPoints

fittingfn = @runICP;
distfn = @evaluateRegistration;
degenfn = @toolinear

% get points from pointclouds
if ~exist('randsample', 'file')
    ind1 = randomsample(cloudsize, sz1);
    ind2 = randomsample(cloudsize, sz2);
else
    ind1 = randsample(cloudsize, sz1);
    ind2 = randsample(cloudsize, sz2);
end

[Model, inliers] = ransacUnknownCorrespondence([pointCloud1(:,ind1);pointCloud2(:,ind2)], fittingfn, distfn, degenfn, minPoints, inlierRadius)


end

function [R,T] = eq_point(x)

[rows, npoints] = size(x);
q = x(1:rows/2,:);
p = x(rows/2+1,:);

m = size(p,2);
n = size(q,2);

% normalize weights
weights = weights ./ sum(weights);

% find data centroid and deviations from centroid
q_bar = q * transpose(weights);
q_mark = q - repmat(q_bar, 1, n);
% Apply weights
q_mark = q_mark .* repmat(weights, 3, 1);

% find data centroid and deviations from centroid
p_bar = p * transpose(weights);
p_mark = p - repmat(p_bar, 1, m);
% Apply weights
%p_mark = p_mark .* repmat(weights, 3, 1);

N = p_mark*transpose(q_mark); % taking points of q in matched order

[U,~,V] = svd(N); % singular value decomposition

R = V*diag([1 1 det(U*V')])*transpose(U);

T = q_bar - R*p_bar;

function r = toolinear(x)
r = true;
return;
% test for planarity
%[linModel1, dist1] = fitline(pointCloud1(1:2,:));
%[linModel2, dist2] = fitline(pointCloud2(1:2,:));
%if (std(dist1) < planeFitParam || std(dist2) < planeFitParam )
%    fprintf('degenerate geometry: too planar for match')
%    % possibly can still do heading matches
%    return;
%end

end
