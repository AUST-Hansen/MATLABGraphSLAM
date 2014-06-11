function PDout = extractSonarFeatures(PD)

%% Gather statistics


% Scanwise curvatures
skip = 50;
[ir,jr] = find(PD.c_i_t);
CurvatureBins = 0:.01:1.1;
ScanWiseCurvatures = zeros(1,length(CurvatureBins));

for iScan = 1:skip:length(PD.timeSteps)
    % Extract scan
    Scan = PD.rangeMeasurements(:,PD.c_i_t(PD.c_i_t(:,iScan)~=0,iScan));
    VehiclePitch = atan(mean(Scan(1,:)./Scan(3,:)));
    % rotate by vehicle pitch (undoing MBsystem stuff)
    ScanRot = Euler2RotMat(0,-VehiclePitch,0)*Scan;
%     ydot = gradient(ScanRot(2,:));
%     yddot = gradient(ydot);
%     zdot = gradient(ScanRot(3,:));
%     zddot = gradient(zdot);
%     kappa = abs(ydot.*zddot - zdot.*yddot)./(ydot.^2 + zdot.^2).^(1.5);    

    zdot = gradient(ScanRot(3,:),ScanRot(2,:));
    zddot = gradient(zdot);
    
    kappa = abs(zddot)./(1+zdot.^2).^(1.5);
    
%     if(~isempty(Scan))
%         subplot(2,1,1)
%         scatter3(ScanRot(2,:),ScanRot(1,:),-ScanRot(3,:))
%         %axis equal
%         axis([-65 -10 -30 30])
%         view(0,0)
%         subplot(2,1,2)
%         %scatter3(PD.Pos(2,:),PD.Pos(1,:),-PD.Pos(3,:),ones(1,size(PD.Pos(2,:),2)));
%         histo = hist(kappa,CurvatureBins);
%         histo = histo./sum(histo);
%         Y = 45*exp(-15*CurvatureBins);
%         ScanWiseCurvatures = ScanWiseCurvatures + histo;
%             plot(CurvatureBins,ScanWiseCurvatures)
%             hold on
%             plot(CurvatureBins,Y,'r')
%         hold off
%         drawnow()
%         pause(.01)
%     end %if

end
% Hard-code thresholds for now based on above analysis
MinCurve = .4;
MaxCurve = .8;
    



% Patch curvatures
for ii = 1:length(PD.timeSteps);
    
        cloud = buildMapFromProcessedData(PD);
        [normals, curvatures] = extractCloudNormals(cloud);
        histo = hist(curvatures,CurvatureBins);
        histo = histo./sum(histo);
        figure; plot(CurvatureBins,histo)
        keyboard

end % for

PDout = [];
