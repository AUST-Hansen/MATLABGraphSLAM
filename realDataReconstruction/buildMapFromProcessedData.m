function pointCloud = buildMapFromProcessedData(PD,tStart, tEnd)


pointCloud = zeros(size(PD.rangeMeasurements));
pointRef = 0;

for ii = tStart:tEnd
    
    R = Euler2RotMat(PD.Phi(ii),PD.Theta(ii),PD.Psi(ii));
    scan_t = PD.rangeMeasurements(:,PD.c_i_t(PD.c_i_t(:,ii)~=0,ii));
    
    pointCloud(:,pointRef+1:pointRef+size(scan_t,2)) = repmat(PD.Pos(:,ii),1,size(scan_t,2)) + R*scan_t;
    pointRef = pointRef+size(scan_t,2);
    
end

pointCloud = pointCloud(:,1:pointRef-1);
    



