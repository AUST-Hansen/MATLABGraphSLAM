function [Rbest, Tbest, inliers_best] = ransacPointCloud(points1,points2,matches)

% matches s.t. p1(:,i) = p2(:,matches(i));
% inliers s.t. p1(:,i) = p2(:,inliers(i));

points1 = points1(:,matches>0)
matches = matches(matches>0);

N = size(points1,2);
population = 1:N
nSamp = int16(5);

maxAddInliers = 0;
maxIter = 10000;
inlierThreshold = 2;
Rbest = [];
Tbest = [];

for t = 1:maxIter
    
    % randomly sample matches
    y = randsample(population,nSamp);
    % and extract
    p1 = points1(:,y);
    p1other = points1(:,~y);
    inliers = zeros(size(matches));
    % left out data for ransacking
    p2 = points2(:,matches(y));
    p2other = points2(:,matches(~y));
    % find data centroid and deviations from centroid
    p1o = mean(p1,2);
    p1_mark = p1 - repmat(p1o, 1, nSamp);
    p2o = mean(p2,2);
    p2_mark = p2 - repmat(p2o, 1, nSamp);

    % get model
    N = p2_mark*transpose(p1_mark); % taking points of q in matched order

    [U,~,V] = svd(N); % singular value decomposition
    
    R = V*diag([1 1 det(U*V')])*transpose(U);
    
    T = p1o - R*p2o;
    
    % test against left out data
    numinliers = 0;
    pPredicted = R*p2 + repmat(T,1,size(p2,2));
    delta = p1 - pPredicted;
    for ii = 1:size(delta,2)
       if norm(delta(:,ii)) < inlierThreshold
           % declare inlier 
           inliers(ii) = 1;
       end
    end
    addInliers = sum(inliers) 
    if addInliers > maxAddInliers
       maxAddInliers = addInliers;
       Rbest = R;
       Tbest = T;
       inliers_best = matches.*inliers;
    end
    
end