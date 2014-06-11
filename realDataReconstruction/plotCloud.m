function plotCloud(pointCloud,fighandle,color)

figure(fighandle);
hold on
        whitebg(fighandle,'k')
        grid off
        scatter3(pointCloud(2,:),pointCloud(1,:),-pointCloud(3,:),2*ones(1,length(pointCloud)),color)
        
        axis equal
        view(-165,6)
