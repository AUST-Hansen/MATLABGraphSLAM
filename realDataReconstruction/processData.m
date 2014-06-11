function ProcessedData = processData(ExtractedData,tStart,tEnd,subMapHalfWidth)
    iStart = find(ExtractedData.timeSteps >= tStart,1);
    iEnd = find(ExtractedData.timeSteps <= tEnd,1,'last');
    
    % reduce timestep
    ProcessedData.timeSteps = ExtractedData.timeSteps(iStart:iEnd) - ExtractedData.timeSteps(iStart);
    % 
    ProcessedData.Pos = ExtractedData.Pos(:,iStart:iEnd);
    ProcessedData.Psi = ExtractedData.Psi(iStart:iEnd);
    ProcessedData.Theta = ExtractedData.Theta(iStart:iEnd);
    ProcessedData.Phi = ExtractedData.Phi(iStart:iEnd);
    ProcessedData.Speed = ExtractedData.Speed(iStart:iEnd);
    ProcessedData.DVL = 0*ProcessedData.Pos;
    %% For now, spoof DVL data
    Xdot = (diff(ProcessedData.Pos')')./repmat(diff(ProcessedData.timeSteps'),3,1);
    for ii = 1:size(ProcessedData.timeSteps,1)-1
        R1_i_v = Euler2RotMat(ProcessedData.Phi(ii),ProcessedData.Theta(ii),ProcessedData.Psi(ii));
        R2_i_v = Euler2RotMat(ProcessedData.Phi(ii+1),ProcessedData.Theta(ii+1),ProcessedData.Psi(ii+1));
        % Average over attitudes at t, t+1 to get body velocity
        u1 = R1_i_v'*Xdot(:,ii);
        u2 = R2_i_v'*Xdot(:,ii);
        ProcessedData.DVL(:,ii) = .5*(u1+u2);
    end
    ProcessedData.DVL(:,end) = ProcessedData.DVL(:,end-1);
    ProcessedData.DVL(1,:) = smooth(ProcessedData.DVL(1,:),10);
    ProcessedData.DVL(2,:) = smooth(ProcessedData.DVL(2,:),10);
    %%
    % Build cleaned cloud
    c_i_t_clean = zeros(512,length(ProcessedData.timeSteps));
    ProcessedData.rangeMeasurements = [];
    rangeCount = 0;
    for ii = iStart:iEnd
        % clean scan
        scan = ExtractedData.rangeMeasurements(:,ExtractedData.c_i_t(ExtractedData.c_i_t(:,ii)~=0,ii));
        rejectEdges = 80;
        cleanedscan = cleanScan(scan(:,rejectEdges:end-rejectEdges),50);
        % record cleaned scan
        ProcessedData.rangeMeasurements = [ProcessedData.rangeMeasurements cleanedscan];   % SLOW!
        c_i_t_clean(1:size(cleanedscan,2),ii-iStart+1) = (rangeCount+1:rangeCount+size(cleanedscan,2))';
        rangeCount = rangeCount + size(cleanedscan,2);
        
    end
    
    % Submap handling
    subMapRefFlag = zeros(1,size(ProcessedData.timeSteps,2));
    subMapIndexRange = zeros(2,size(ProcessedData.timeSteps,2));
    submapIndex = 1;

    for ii = 1 + subMapHalfWidth:subMapHalfWidth:size(ProcessedData.timeSteps,1)-subMapHalfWidth
       subMapRefFlag(ii) = 1;
       subMapIndexRange(1,ii) = c_i_t_clean(1,ii-subMapHalfWidth);
       subMapIndexRange(2,ii) = max(max(c_i_t_clean(:,ii:ii+subMapHalfWidth)));
       submapIndex = submapIndex + 1; 
    end
    
    ProcessedData.c_i_t = c_i_t_clean;
    ProcessedData.subMapRefFlag = subMapRefFlag;
    ProcessedData.subMapIndexRange = subMapIndexRange;