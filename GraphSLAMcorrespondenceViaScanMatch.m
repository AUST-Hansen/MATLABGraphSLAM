function [correspondences_new,relHeading] = GraphSLAMcorrespondenceViaScanMatch(mu,correspondences,posThreshold,probThreshold,RMStolerance,measIndices,rangeMeasurements,imagenex,imagenexSkip,LCWeight,InlineWeight,varargin)

DEBUG = false;
stateHist = [];

if(~isempty(varargin))
    for k = 1:nVarargs
        switch(varargin{k})
            case{'State'}
                stateHist = varargin{k+1};
            case{'DEBUG'}
                DEBUG = true;
        end %switch
    end %for
end %if

% initialize correspondences.c_i_t output
correspondences_new = correspondences;
Nstates = size(correspondences.c_i_t,2);
stateSize = 6; % bad practice
size(mu(stateSize*Nstates+1:end),1)/3;
stateHist = reshape(mu(1:stateSize*Nstates),stateSize,[]);
mapFeatures = reshape(mu(stateSize*Nstates+1:end),3,[]);
relHeading.validMatch = spalloc(Nstates,Nstates,6*Nstates);
relHeading.deltaHeading = spalloc(Nstates,Nstates,6*Nstates);

endLCidx = 0; % last index you'll allow loop closure scan matching.

ignxEndIdx = floor(imagenex.numBeams/imagenexSkip);

for ii = 1:Nstates
    
    pos1 = stateHist(1:2,ii);
    heading1 = stateHist(5,ii);
    berg_R_veh1 = [cos(heading1), -sin(heading1) 0; sin(heading1), cos(heading1), 0;0,0,1];
    pose1FeatureIndex = find (correspondences.c_i_t(:,ii)~=-17);
    featuresSeenAtPose1 = correspondences.c_i_t(pose1FeatureIndex,ii);
    %pointCloud1 = mapFeatures(:,featuresSeenAtPose1);
    pointCloud1 = rangeMeasurements(:,measIndices(measIndices(1:ignxEndIdx,ii)~=-17,ii));
    pointCloudA = berg_R_veh1*pointCloud1 + repmat([pos1;0],1,size(pointCloud1,2));
    pointCloud1 = pointCloudA;
    for jj = ii+1:Nstates
        goodmatch = false;
        if (sum(correspondences.c_i_t(1:ignxEndIdx,ii) ~= -17) == 0)
            continue
        end
        pos2 = stateHist(1:2,jj);
        heading2 = stateHist(5,jj);
        berg_R_veh2 = [cos(heading2), -sin(heading2) 0; sin(heading2), cos(heading2), 0;0,0,1];
        % are they close-ish?
        % if they're close and EITHER not far apart in time, or we're
        % looking for loop closure with first few soundings
        if (norm(pos1-pos2) < posThreshold )
            if (sum(correspondences.c_i_t(:,jj) ~= -17) == 0)
                continue
            end
            pose2FeatureIndex = find (correspondences.c_i_t(:,jj)~=-17);
            featuresSeenAtPose2 = correspondences.c_i_t(pose2FeatureIndex,jj);
            %pointCloud2 = mapFeatures(:,featuresSeenAtPose2);
            pointCloud2 = rangeMeasurements(:,measIndices(measIndices(1:ignxEndIdx,jj)~=-17,jj));
            pointCloudB = berg_R_veh2*pointCloud2 + repmat([pos2;0],1,size(pointCloud2,2));
            pointCloud2 = pointCloudB;
            % Do ICP
            %[R,T,ERR] = icp(pointCloud1,pointCloud2,20,'WorstRejection',.1,'twoDee',true);
            [R,T,ERR,idxKNN,dKNN] = robustScanMatch(pointCloud1,pointCloud2);
            %ERR;
            pNew = R*pointCloud2 + repmat(T,1,size(pointCloud2,2));
            % Do KNN
            %[idxKNN, dKNN]=knnsearch(pointCloud1',pNew');
            
            ERR(end)
            
            [ii jj]
            for kk = 1:length(idxKNN)
                if (~isempty(ERR) && dKNN(kk) <= probThreshold && ERR(end) < RMStolerance)
                    %fprintf('rms: %d\n',ERR(end));
                    % looks like we have a winner
                    % kk is the index in pointCloud2
                    % idxKNN(kk) is the corresponding index in pointCloud1
                    %
                    if (isempty(   intersect(correspondences_new.c_i_t(:,jj), correspondences_new.c_i_t(pose1FeatureIndex(idxKNN(kk)),ii)))) % mutex constraint
                        correspondences_new.c_i_t(pose2FeatureIndex(kk),jj) = correspondences_new.c_i_t(pose1FeatureIndex(idxKNN(kk)),ii);
                        if (abs(jj - ii) > 600)
                            correspondences_new.weight(pose2FeatureIndex(kk),jj) = LCWeight;
                        else
                            correspondences_new.weight(pose2FeatureIndex(kk),jj) = InlineWeight;
                        end
                        %fprintf('match between %d and %d\n',ii,jj)
                        goodmatch = true;
                    end
                end
            end
            if (goodmatch)
                relHeading.validMatch(ii,jj) = 1;
                relHeadingMat = berg_R_veh1'*R*berg_R_veh2;
                relHeading.deltaHeading(ii,jj) = -asin(relHeadingMat(1,2));
            end
        end
        
        %%%%%% DEBUG stuff
        if (DEBUG && rand()< 1  && goodmatch == true  )
            figure(9)
            hold on; axis equal;
            for qq = 1:length(featuresSeenAtPose1)
                scatter(pointCloud1(1,qq),pointCloud1(2,qq),'r')
                %text(pointCloud1(1,qq),pointCloud1(2,qq)+1,num2str(featuresSeenAtPose1(qq)));
            end
            for qq= 1:length(featuresSeenAtPose2)
                scatter(pNew(1,qq),pNew(2,qq),'b')
                %text(pNew(1,qq),pNew(2,qq)-1,num2str(featuresSeenAtPose2(qq)));
            end
            keyboard
            close(9)
        end
        
        
    end
    

    
    
    
    
end
    % Retire unseen features
UniqueFeatures = unique(correspondences_new.c_i_t(correspondences_new.c_i_t>0))
for iNew = 1:length(UniqueFeatures)
    correspondences_new.c_i_t(correspondences_new.c_i_t==UniqueFeatures(iNew)) = iNew;
end
