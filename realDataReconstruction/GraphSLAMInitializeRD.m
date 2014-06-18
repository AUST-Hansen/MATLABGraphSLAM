function [SLAMdata, OtherData] = GraphSLAMInitializeRD(PD)

DEBUG = true;
corruptDataWithMotion = false;
xy = 0*PD.Pos(1:2,:);
uv = PD.DVL(1:2,:);
for ii = 2:size(PD.timeSteps,1);
    
end
% Deal with wrapping
psi0 = PD.Psi';
wrap = 0;
psi(1) = psi0(1);
for ii = 2:length(psi0);
    if(psi0(ii-1) > 6 && psi0(ii) < 1 ) % positive wrap has occurred
        wrap = wrap+1;
    elseif (psi0(ii-1) < 1 && psi0(ii) >6 ) % negative wrap has occurred
        wrap = wrap-1;
    end
    psi(ii) = psi0(ii) + 2*pi*wrap;
end
psi = psi - psi(1);

omega = gradient(psi,PD.timeSteps(2)-PD.timeSteps(1));
psiRedo = 0*psi;
for jj = 2:length(psiRedo)
    psiRedo(jj) = psiRedo(jj-1) + omega(jj-1)*(PD.timeSteps(jj) - PD.timeSteps(jj-1));
end
figure; plot(psi); hold on; plot(omega,'g'),plot(psiRedo,'r')

bias = 0*psi;

if (corruptDataWithMotion) % adding apparent iceberg drift
    psiBerg = .00012*PD.timeSteps;
    for ii = 2:size(PD.timeSteps,1);
        R = Euler2RotMat(0,0,psi(ii-1) - psiBerg(ii-1));
        xy(:,ii) = xy(:,ii-1) + R(1:2,1:2)*uv(1:2,ii-1)*(PD.timeSteps(ii)-PD.timeSteps(ii-1));
    end
else
    for ii = 2:size(PD.timeSteps,1);
        R = Euler2RotMat(0,0,psi(ii-1));
        xy(:,ii) = xy(:,ii-1) + R(1:2,1:2)*uv(1:2,ii-1)*(PD.timeSteps(ii)-PD.timeSteps(ii-1));
    end
end
% Initialize slam trajectory quantities
trajectory = reshape([xy;uv;psi;bias],[],1);
SLAMdata.c_i_t = PD.c_i_t;
SLAMdata.stateSize = 6;
SLAMdata.Nstates = size(PD.timeSteps,1);
SLAMdata.mapDimension = 4; % xyz + z_normal
% Initialize map features
PDnew = PD;
PDnew.Pos(1:2,:) = xy;
PDnew.Psi = psi';
pointCloud = buildFeatureMapFromProcessedData(PDnew,1,size(PD.timeSteps,1));
size(pointCloud)
size(PDnew.Descriptors(4,:))
pointCloudWithZNormal = [pointCloud;PDnew.Descriptors(4,:)];
SLAMdata.mu = [trajectory; reshape(pointCloudWithZNormal,[],1)];
% Other quantities
OtherData.Descriptors = PD.Descriptors;
    OtherData.Z = PD.Pos(3,:);
    OtherData.timeSteps = PD.timeSteps;
    OtherData.inputs = [1.5*ones(size(PD.timeSteps'));omega];
    OtherData.DVL = PD.DVL;
    OtherData.Theta = PD.Theta;
    OtherData.Phi = PD.Phi;
    OtherData.measIndex = PD.c_i_t;
   
    if(DEBUG)
        figure;
        statehist = reshape(SLAMdata.mu(1:SLAMdata.stateSize*SLAMdata.Nstates),6,[]);
        scatter3(statehist(2,:),statehist(1,:),-OtherData.Z)
        axis equal
        hold on
        points = reshape(SLAMdata.mu(SLAMdata.stateSize*15570+1:end),SLAMdata.mapDimension,[]);
        scatter3(points(2,:),points(1,:),-points(3,:),3*ones(1,size(points,2)),'g')
        view(130,16)
        title('Initial estimate')
    end