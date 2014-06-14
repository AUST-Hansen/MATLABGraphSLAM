function PDout = sweepSonarFeatures(PD)

   subMapHalfWidth = 20;
   tStart = 1;
   tEnd = length(PD.timeSteps);
   tSkip = 1;
   PLOTS = false;
   % Max and min curvatures: hard-coded
   MinCurve = .6;
   MaxCurve = 1.;
   min2dCurve = .03; 
   max2dCurve = .05;
   % Build initial submap
   subMap = [];
   %% feature extraction bins
   c_i_t_features = zeros(20,tEnd-tStart+1);
   descriptors = [];
   featureCount = 0;
   while(isempty(subMap) || size(subMap,2) < 5000)
    subMap = buildMapFromProcessedData(PD,tStart,tStart+2*subMapHalfWidth+1);
    tStart = tStart + 1
   end
   
   
   size(subMap)
   
   tStart = tStart - 1;
   if (PLOTS)
       figure(5)
       view(-165,6)
       whitebg(5,'k')
   end
   for iScan = subMapHalfWidth+tStart:tSkip:tEnd - subMapHalfWidth
       subMap = buildMapFromProcessedData(PD,iScan-subMapHalfWidth,iScan+subMapHalfWidth);
       if (size(subMap,2) < 5000)
           continue
       end
       % Extract scan
       Scan = PD.rangeMeasurements(:,PD.c_i_t(PD.c_i_t(:,iScan)~=0,iScan));
       Riv_ref = Euler2RotMat(PD.Phi(iScan),PD.Theta(iScan),PD.Psi(iScan));
       % World frame
       ScanWF = repmat(PD.Pos(:,iScan),1,size(Scan,2)) + Riv_ref*Scan;

    if(PLOTS)
       plotCloud(subMap);
       hold on;
       scatter3(ScanWF(2,:),ScanWF(1,:),-ScanWF(3,:),5*ones(1,length(ScanWF)),'b+')
    end
       VehiclePitch = atan(mean(Scan(1,:)./Scan(3,:)));
       % rotate by vehicle pitch (undoing MBsystem stuff)
       ScanRot = Euler2RotMat(0,-VehiclePitch,0)*Scan;
       zdot = gradient(ScanRot(3,:),ScanRot(2,:));
       zddot = gradient(smooth(zdot)');
       % curvature is kappa
       kappa = abs(zddot)./(1+zdot.^2).^(1.5);
       % Extract points with high curvature in scanwise direction for further
       % processing
       CurvyPoints = ScanWF(:,kappa>MinCurve & kappa<MaxCurve);
       CurvyPointsVF = Scan(:,kappa>MinCurve & kappa<MaxCurve);
       kappaCandidates = kappa(kappa>MinCurve & kappa<MaxCurve);
       if(~isempty(CurvyPoints) && PLOTS)
           scatter3(CurvyPoints(2,:),CurvyPoints(1,:),-CurvyPoints(3,:),100*ones(1,size(CurvyPoints,2)),'g.')
       end
       
       %% Now do more processing on curvy points
       neighborIndices = rangesearch(subMap',CurvyPoints',5);
       curvatures = zeros(1,size(CurvyPoints,2));
       normals = zeros(3,size(CurvyPoints,2));
       for ii = 1:size(CurvyPoints,2)
           
           % weed out borders and spurious points
           if(length(neighborIndices{ii}) < 10)
               continue;
           end
           
           % make neighbor matrix
           P = subMap(:,neighborIndices{ii});
           Po = mean(P,2);
           Preg = P - repmat(Po,1,size(P,2));
           M = Preg*Preg';
           [u s v] = svd(M);
           normal = v(:,3);
           if (normal'*(Po-PD.Pos(:,iScan)) > 0)
               normal = -normal;
           end
           normal = Riv_ref'*normal;
           normals(:,ii) = normal;
           curvatures(ii) = s(3,3)/(s(3,3)+s(2,2)+s(1,1)); %[0 0 1]*v(:,3)*sign([1 0 0]*v(:,3));
       end
       
       % filter out low curvatures

       CurvyPoints2D = CurvyPoints(:,curvatures>min2dCurve&curvatures<max2dCurve);
       CurvyPoints2DVF = CurvyPointsVF(:,curvatures>min2dCurve&curvatures<max2dCurve);
       kappaCulled = kappaCandidates(curvatures>min2dCurve&curvatures<max2dCurve);
       curvesCulled = curvatures(curvatures>min2dCurve&curvatures<max2dCurve);
       normalsCulled = normals(:,curvatures>min2dCurve&curvatures<max2dCurve);

       
       if(~isempty(CurvyPoints) && PLOTS)
           figure(1);
           plot(curvatures)
           axis([0 size(CurvyPoints,2) 0 .1])
           figure(5)
           scatter3(CurvyPoints(2,:),CurvyPoints(1,:),-CurvyPoints(3,:),100*ones(1,size(CurvyPoints,2)),'g.')
           scatter3(CurvyPoints2D(2,:),CurvyPoints2D(1,:),-CurvyPoints2D(3,:),1000*ones(1,size(CurvyPoints2D,2)),'c.')
       end
       if(PLOTS)
           drawnow()
           pause(.01)
           hold off;
       end
       
       % Build feature descriptor
       % Position in vehicle frame 
       % nz component
       % normal grazing angle
       normalHeading = mod(atan2(-normalsCulled(1,:),-normalsCulled(2,:)),2*pi);
       grazingAngles = normalHeading - repmat(PD.Psi(iScan),1,size(normalHeading,2)) + pi/2;
       % 1D curvature
       % 2d curvature

       if (size(CurvyPoints2DVF,2) > 0)
           descriptors = [descriptors [CurvyPoints2DVF;...
               normalsCulled(3,:);...
               grazingAngles;...
               kappaCulled;...
               curvesCulled]];
           c_i_t_features(1:size(curvesCulled,2),iScan) = (featureCount+1:featureCount+size(curvesCulled,2))';
           featureCount = featureCount + size(curvesCulled,2);
       end
       [iScan size(curvesCulled,2)]
   end %for iScan
   
   PDout.timeSteps = PD.timeSteps;
   PDout.Pos = PD.Pos;
   PDout.Psi = PD.Psi;
   PDout.Theta = PD.Theta;
   PDout.Phi = PD.Phi;
   PDout.Speed = PD.Speed;
   PDout.DVL = PD.DVL;
   PDout.Descriptors = descriptors;
   PDout.c_i_t = c_i_t_features;