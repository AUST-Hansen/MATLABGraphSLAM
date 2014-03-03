function c_i_t_new = lookForLoopClosure(c_i_t,rangeMeasurements,measurementTimestamps,probThreshold)

jj=1;
offset = 650;
c_i_t_new = c_i_t;
pose1FeatureIndex = find (c_i_t(:,jj)~=-17);
featuresSeenAtPose1 = c_i_t(pose1FeatureIndex,jj);
pointCloud1 = rangeMeasurements(:,featuresSeenAtPose1);
resids = 1000*ones(size(measurementTimestamps));

for ii = jj+offset:size(measurementTimestamps,2)
    % build poses
    pose2FeatureIndex = find (c_i_t(:,ii)~=-17);
    if (isempty(pose2FeatureIndex))
        continue
    end
    featuresSeenAtPose2 = c_i_t(pose2FeatureIndex,ii);
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
pose2FeatureIndex = find (c_i_t(:,ii)~=-17);
featuresSeenAtPose2 = c_i_t(pose2FeatureIndex,ii);
pointCloud2 = rangeMeasurements(:,featuresSeenAtPose2);
% do icp
[R,T,ERR,idxKNN,dKNN] = robustScanMatch(pointCloud1,pointCloud2);
pNew = R*pointCloud2 + repmat(T,1,size(pointCloud2,2));
% Do KNN
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

    % Retire unseen features
UniqueFeatures = unique(c_i_t_new(c_i_t_new>0));
for iNew = 1:length(UniqueFeatures)
    c_i_t_new(c_i_t_new==UniqueFeatures(iNew)) = iNew;
end

%%%%%% DEBUG stuff
if (true)
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
        