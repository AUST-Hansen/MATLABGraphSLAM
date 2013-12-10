function c_i_t_new = lookForLoopClosure(c_i_t,rangeMeasurements,measurementTimestamps,probThreshold)

jj=1;
c_i_t_new = c_i_t;
pose1FeatureIndex = find (c_i_t(:,jj)~=-17);
featuresSeenAtPose1 = c_i_t(pose1FeatureIndex,jj);
pointCloud1 = rangeMeasurements(:,featuresSeenAtPose1);
resids = 0*measurementTimestamps;

for ii = jj+1:size(measurementTimestamps,2)
    % build poses
    pose2FeatureIndex = find (c_i_t(:,ii)~=-17);
    featuresSeenAtPose2 = c_i_t(pose2FeatureIndex,ii);
    pointCloud2 = rangeMeasurements(:,featuresSeenAtPose2);
    % do icp
    [R,T,ERR] = icp(pointCloud1,pointCloud2);
    
    resids(ii) = ERR(end);
    if(mod(ii,50) == 0)
        %fprintf('%d of %d poses\n',ii,size(measurementTimestamps,2)) ;
    end
end

lcResids = resids;
lcResids(1:50) = 1000;
%% look for loop closure
[y,ii]=min(lcResids);
pose2FeatureIndex = find (c_i_t(:,ii)~=-17);
featuresSeenAtPose2 = c_i_t(pose2FeatureIndex,ii);
pointCloud2 = rangeMeasurements(:,featuresSeenAtPose2);
% do icp
[R,T,ERR] = icp(pointCloud1,pointCloud2);
pNew = R*pointCloud2 + repmat(T,1,size(pointCloud2,2));
% Do KNN
[idxKNN, dKNN]=knnsearch(pointCloud1',pNew');
goodmatch = false;
for kk = 1:length(idxKNN)
    if (dKNN(kk) <= probThreshold)
        if (isempty(   intersect(c_i_t_new(:,ii), c_i_t_new(pose1FeatureIndex(idxKNN(kk)),jj)))) % mutex constraint
            c_i_t_new(pose2FeatureIndex(kk),ii) = c_i_t_new(pose1FeatureIndex(idxKNN(kk)),jj);
            fprintf('match between %d and %d\n',jj,ii)
            goodmatch = true;
        end
    end
end