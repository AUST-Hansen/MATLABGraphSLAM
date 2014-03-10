function correspondences_new = lookForLoopClosure(correspondences,rangeMeasurements,measurementTimestamps,probThreshold,LCWeight,varargin)

jj=1;
offset = 680;
endidx = size(measurementTimestamps,2);
correspondences_new = correspondences;
pose1FeatureIndex = find (correspondences.c_i_t(:,jj)~=-17);
featuresSeenAtPose1 = correspondences.c_i_t(pose1FeatureIndex,jj);
pointCloud1 = rangeMeasurements(:,featuresSeenAtPose1);
resids = 1000*ones(size(measurementTimestamps));
DEBUG = false;
% Parse varargin
nVarargs = length(varargin);
if(~isempty(varargin))
    for k = 1:nVarargs
        switch(varargin{k})
            case{'MinIndex'}
                offset = varargin{k+1};
            case{'MaxIndex'}
                endidx = varargin{k+1};

            case{'DEBUG'}
                DEBUG = true;
        end %switch
    end %for
end %if

for ii = jj+offset:endidx
    % build poses
    pose2FeatureIndex = find (correspondences.c_i_t(:,ii)~=-17);
    if (isempty(pose2FeatureIndex))
        continue
    end
    featuresSeenAtPose2 = correspondences.c_i_t(pose2FeatureIndex,ii);
    pointCloud2 = rangeMeasurements(:,featuresSeenAtPose2);
    % do icp
    [R,T,ERR,idxKNN,dKNN] = robustScanMatch(pointCloud1,pointCloud2);
    %[R,T,ERR] = icp(pointCloud1,pointCloud2);  
    resids(ii) = ERR(end);
    if(mod(ii,50) == 0)
        %fprintf('%d of %d poses\n',ii,size(measurementTimestamps,2)) ;
    end
    ii
end

lcResids = resids;
lcResids(1:50) = 1000;
%% look for loop closure
[y,ii]=min(lcResids);
pose2FeatureIndex = find (correspondences.c_i_t(:,ii)~=-17);
featuresSeenAtPose2 = correspondences.c_i_t(pose2FeatureIndex,ii);
pointCloud2 = rangeMeasurements(:,featuresSeenAtPose2);
% do icp
[R,T,ERR,idxKNN,dKNN] = robustScanMatch(pointCloud1,pointCloud2,'MaxIterations',200);
pNew = R*pointCloud2 + repmat(T,1,size(pointCloud2,2));
% Do KNN
goodmatch = false;
for kk = 1:length(idxKNN)
    if (dKNN(kk) <= probThreshold)
        if (isempty(   intersect(correspondences_new.c_i_t(:,ii), correspondences_new.c_i_t(pose1FeatureIndex(idxKNN(kk)),jj)))) % mutex constraint
            correspondences_new.c_i_t(pose2FeatureIndex(kk),ii) = correspondences_new.c_i_t(pose1FeatureIndex(idxKNN(kk)),jj);
            fprintf('match between %d and %d\n',jj,ii)
            goodmatch = true;
            correspondences_new.weight(pose2FeatureIndex(kk),ii) = LCWeight;
        end
    end
end

    % Retire unseen features
UniqueFeatures = unique(correspondences_new.c_i_t(correspondences_new.c_i_t>0));
for iNew = 1:length(UniqueFeatures)
    correspondences_new.c_i_t(correspondences_new.c_i_t==UniqueFeatures(iNew)) = iNew;
end

%%%%%% DEBUG stuff
if (DEBUG)
    figure(9)
    hold on; axis equal;
    for qq = 1:length(featuresSeenAtPose1)
        scatter(pointCloud1(1,qq),pointCloud1(2,qq),'r')
        text(pointCloud1(1,qq),pointCloud1(2,qq)+1,num2str(featuresSeenAtPose1(qq)));
    end
    for qq= 1:length(featuresSeenAtPose2)
        scatter(pNew(1,qq),pNew(2,qq),'b')
        %scatter(pointCloud2(1,qq),pointCloud2(2,qq),'g')
        text(pNew(1,qq),pNew(2,qq)-1,num2str(featuresSeenAtPose2(qq)));
    end
      
    keyboard
    close(9)
end
        