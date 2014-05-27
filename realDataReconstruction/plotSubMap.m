function pointCloud = plotSubMap(PD,index)

if(PD.subMapRefFlag(index) == 0)
    fprintf('no submap associated with index %d\n',index)
    fprintf('available indices:\n')
    ind = find(PD.subMapRefFlag);
    fprintf(num2str(ind));
    fprintf('\n');
    pointCloud = [];
else
    minIndex = PD.subMapIndexRange(1,index);
    maxIndex = PD.subMapIndexRange(2,index);
    [~,jmin] = find(PD.c_i_t == minIndex);
    [~,jmax] = find(PD.c_i_t == maxIndex);
    pointCloud = zeros(3,maxIndex-minIndex+1);
    pointRef = 0;
    % build map
    for ii = jmin:jmax
        ii
       R = Euler2RotMat(PD.Phi(ii),PD.Theta(ii),PD.Psi(ii));
       scan_t = PD.rangeMeasurements(:,PD.c_i_t(PD.c_i_t(:,ii)~=0,ii));
       
       pointCloud(:,pointRef+1:pointRef+size(scan_t,2)) = repmat(PD.Pos(:,ii),1,size(scan_t,2)) + R*scan_t; 
       pointRef = pointRef+size(scan_t,2);
        
    end
    
    
    scatter3(pointCloud(2,:),pointCloud(1,:),-pointCloud(3,:),1*ones(1,length(pointCloud)),'b')
    
    hold on
    scatter3(PD.Pos(2,index),PD.Pos(1,index),-PD.Pos(3,index),'g')
    axis equal
    view(-148,6)
    hold off
end

