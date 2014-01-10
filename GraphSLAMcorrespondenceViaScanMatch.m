function [c_i_t_new] = GraphSLAMcorrespondenceViaScanMatch(mu,c_i_t,posThreshold,probThreshold,RMStolerance,measIndices,rangeMeasurements,imagenex,imagenexSkip)

% initialize c_i_t output
c_i_t_new = c_i_t;
Nstates = size(c_i_t,2);
stateSize = 6; % bad practice
size(mu(stateSize*Nstates+1:end),1)/3;
stateHist = reshape(mu(1:stateSize*Nstates),stateSize,[]);
mapFeatures = reshape(mu(stateSize*Nstates+1:end),3,[]);

for ii = 1:Nstates
    pos1 = stateHist(1:2,ii);
    heading1 = stateHist(5,ii);
    berg_R_veh1 = [cos(heading1), -sin(heading1) 0; sin(heading1), cos(heading1), 0;0,0,1];
    pose1FeatureIndex = find (c_i_t(:,ii)~=-17);
    featuresSeenAtPose1 = c_i_t(pose1FeatureIndex,ii);
    %pointCloud1 = mapFeatures(:,featuresSeenAtPose1);
    pointCloud1 = rangeMeasurements(:,measIndices(measIndices(1:floor(imagenex.numBeams/imagenexSkip),ii)~=-17,ii));
    pointCloudA = berg_R_veh1*pointCloud1 + repmat([pos1;0],1,size(pointCloud1,2));
    pointCloud1 = pointCloudA;
    for jj = ii+1:Nstates
        if (sum(c_i_t(1:floor(imagenex.numBeams/imagenexSkip),ii) ~= -17) == 0)
            continue
        end
        pos2 = stateHist(1:2,jj);
        heading2 = stateHist(5,jj);
        berg_R_veh2 = [cos(heading2), -sin(heading2) 0; sin(heading2), cos(heading2), 0;0,0,1];
        % are they close-ish?
        % if they're close and EITHER not far apart in time, or we're
        % looking for loop closure with first few soundings
        if ((norm(pos1-pos2) < posThreshold && abs(ii-jj)<30 ) || (ii<20 && norm(pos1-pos2) < posThreshold)  )
            if (sum(c_i_t(:,jj) ~= -17) == 0)
                continue
            end
            pose2FeatureIndex = find (c_i_t(:,jj)~=-17);
            featuresSeenAtPose2 = c_i_t(pose2FeatureIndex,jj);
            %pointCloud2 = mapFeatures(:,featuresSeenAtPose2);
            pointCloud2 = rangeMeasurements(:,measIndices(measIndices(1:floor(imagenex.numBeams/imagenexSkip),jj)~=-17,jj));
            pointCloudB = berg_R_veh2*pointCloud2 + repmat([pos2;0],1,size(pointCloud2,2));
            pointCloud2 = pointCloudB;
            % Do ICP
            [R,T,ERR] = icp(pointCloud1,pointCloud2,'WorstRejection',.1);
            ERR
            pNew = R*pointCloud2 + repmat(T,1,size(pointCloud2,2));
            % Do KNN
            [idxKNN, dKNN]=knnsearch(pointCloud1',pNew');
            goodmatch = false;
            ERR(end)
            for kk = 1:length(idxKNN)
                if (ERR(end) < RMStolerance && dKNN(kk) <= probThreshold)
                    %fprintf('rms: %d\n',ERR(end));
                    % looks like we have a winner
                    % kk is the index in pointCloud2
                    % idxKNN(kk) is the corresponding index in pointCloud1
                    %
                    if (isempty(   intersect(c_i_t_new(:,jj), c_i_t_new(pose1FeatureIndex(idxKNN(kk)),ii)))) % mutex constraint
                        c_i_t_new(pose2FeatureIndex(kk),jj) = c_i_t_new(pose1FeatureIndex(idxKNN(kk)),ii);
                        fprintf('match between %d and %d\n',ii,jj)
                        goodmatch = true;
                    end
                end
            end
        end
        
        %%%%%% DEBUG stuff
        if (false && rand()<.005  && goodmatch == true )
            figure(9)
            hold on; axis equal;
            for qq = 1:length(featuresSeenAtPose1)
                scatter(pointCloud1(1,qq),pointCloud1(2,qq),'r')
                text(pointCloud1(1,qq),pointCloud1(2,qq)+1,num2str(featuresSeenAtPose1(qq)));
            end
            for qq= 1:length(featuresSeenAtPose2)
                scatter(pNew(1,qq),pNew(2,qq),'b')
                text(pNew(1,qq),pNew(2,qq)-1,num2str(featuresSeenAtPose2(qq)));
            end
            keyboard
            close(9)
        end
        
        
    end
    

    
    
    
    
end
    % Retire unseen features
UniqueFeatures = unique(c_i_t_new(c_i_t_new>0))
for iNew = 1:length(UniqueFeatures)
    c_i_t_new(c_i_t_new==UniqueFeatures(iNew)) = iNew;
end
