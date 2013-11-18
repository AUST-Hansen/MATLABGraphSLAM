function [c_i_t_new] = GraphSLAMcorrespondenceViaScanMatch(mu,c_i_t,posThreshold,probThreshold,RMStolerance)

% initialize c_i_t output
c_i_t_new = c_i_t;
Nstates = size(c_i_t,2);
stateSize = 6; % bad practice
size(mu(stateSize*Nstates+1:end),1)/3;
stateHist = reshape(mu(1:stateSize*Nstates),stateSize,[]);
mapFeatures = reshape(mu(stateSize*Nstates+1:end),3,[]);

for ii = 1:Nstates
    pos1 = stateHist(1:2,ii);
    pose1FeatureIndex = find (c_i_t(:,ii)~=-17);
    featuresSeenAtPose1 = c_i_t(pose1FeatureIndex,ii);
    pointCloud1 = mapFeatures(:,featuresSeenAtPose1);
    for jj = ii+1:Nstates
        if (sum(c_i_t(:,ii) ~= -17) == 0)
            continue
        end
        pos2 = stateHist(1:2,jj);
        % are they close-ish?
        if (norm(pos1-pos2) < posThreshold)
            if (sum(c_i_t(:,jj) ~= -17) == 0)
                continue
            end
            pose2FeatureIndex = find (c_i_t(:,jj)~=-17);
            featuresSeenAtPose2 = c_i_t(pose2FeatureIndex,jj);
            pointCloud2 = mapFeatures(:,featuresSeenAtPose2);
            % Do ICP
            [R,T,ERR] = icp(pointCloud1,pointCloud2,'twoDee',true);
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
                    %if (isempty(   intersect(c_i_t_new(:,jj), c_i_t_new(pose1FeatureIndex(idxKNN(kk)),ii)))) % mutex constraint
                    c_i_t_new(pose2FeatureIndex(kk),jj) = c_i_t_new(pose1FeatureIndex(idxKNN(kk)),ii);
                    fprintf('match between %d and %d\n',ii,jj)
                    goodmatch = true;
                    %end
                end
            end
        end
        
        %%%%%% DEBUG stuff
        if (false  && goodmatch == true)
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