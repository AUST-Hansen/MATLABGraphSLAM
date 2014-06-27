% tuning script
clear all; close all; clc;
%filenames = ['ProcessedData0_2000.mat' 'ProcessedData3500_5200.mat']
load ProcessedData0_2000
if(~exist('PDfeatures0_2000.mat','file'))
    
    PDfeatures = sweepSonarFeatures(ProcessedData);
    save('PDfeatures0_2000.mat','PDfeatures')
else
    load PDfeatures0_2000.mat
end
%%
cloud = buildFeatureMapFromProcessedData(PDfeatures,1,size(PDfeatures.c_i_t,2));
cloud2 = buildMapFromProcessedData(ProcessedData,1,size(ProcessedData.c_i_t,2));
%%
scale = [min(cloud'); max(cloud')];
scale = reshape(scale,1,[]);
scale = [scale(3:4) scale(1:2) -scale(6:5)]
figure(5)
axis manual

hold on;
figure(5)
axis(scale)
subplot(4,1,2)

plotCloud(cloud,5,'r',2)
%%
%figure(3); plot(PDfeatures.Descriptors(4:end,:)'),legend('nz','grazing angle','k1','k2')
PDfeatures1 = PDfeatures;
clear ProcessedData
clear PDfeatures
%% second swath
load ProcessedData3500_5200
if(~exist('PDfeatures3500_5200.mat','file'))
    load ProcessedData3500_5200
    PDfeatures = sweepSonarFeatures(ProcessedData);
    save('PDfeatures3500_5200.mat','PDfeatures')
else
    load PDfeatures3500_5200.mat
end
CORRUPT_DATA = false;
if (CORRUPT_DATA)
    PDfeatures.Pos(1,:) = PDfeatures.Pos(1,:) + 40;
    PDfeatures.Pos(2,:) = PDfeatures.Pos(2,:) + 40;
    PDfeatures.Psi = PDfeatures.Psi + linspace(1,1,size(PDfeatures.Psi,1))';
end
cloud3 = buildFeatureMapFromProcessedData(PDfeatures,1,size(PDfeatures.c_i_t,2));
cloud4 = buildMapFromProcessedData(ProcessedData,1,size(ProcessedData.c_i_t,2));
%%
subplot(4,1,1)

plotCloud(cloud2(:,1:100:end),5,'r',2);
hold on
plotCloud(cloud4(:,1:100:end),5,'b',2);
axis manual
axis(scale)
title('raw data')
hold on;
%%
subplot(4,1,2)
plotCloud(cloud3,5,'b',2)
axis manual
axis(scale)
title('features')
PDfeatures2 = PDfeatures;
clear PDfeatures

% just overlap

% all points in second cloud within five meters of first
ball = 30;
W = diag([.5 .5 20 30]);
hypercloud = W*[cloud; PDfeatures1.Descriptors(4,:)];
hypercloud3 = W*[cloud3; PDfeatures2.Descriptors(4,:)];
% Malanhobis distance metric
[idx dist] = knnsearch(hypercloud3',hypercloud');
 idxInsideRadius = idx(dist<ball);
 matchcloud1 = cloud(:,dist<ball);
 matchcloud2 = cloud3(:,idxInsideRadius);
[idx2 dist2] = knnsearch(hypercloud(:,dist<ball)',hypercloud3(:,idxInsideRadius)');
idxInsideRadius2 = idx2(dist2<ball);
matchcloud1p = matchcloud1(:,idxInsideRadius2);
matchcloud2p = matchcloud2(:,dist2<ball);
% [idx dist] = knnsearch(cloud3',cloud');
% idxCloud1 = 1:length(idx);
% idxCloud1 = idxCloud1(dist<ball);
% idxInsideRadius = idx(dist<ball);
% matchcloud1 = cloud(:,dist<ball);
% matchcloud2 = cloud3(:,idxInsideRadius);
% [idx2 dist2] = knnsearch(matchcloud1',matchcloud2');
% idxCloud2 = 1:length(idx2);
% idxCloud2 = idxCloud2(dist2<ball);
% idxInsideRadius2 = idx2(dist2<ball);
% matchcloud1prime = matchcloud1(:,idxInsideRadius2);
% matchcloud2prime = matchcloud2(:,dist2<ball);
    
normalAgreement = abs(PDfeatures1.Descriptors(4,dist<ball)'-PDfeatures2.Descriptors(4,idxInsideRadius)');
subplot(4,1,3)
plotCloud(matchcloud1,5,'r',5)
hold on
plotCloud(matchcloud2,5,'b',5)
%plotCloud(matchcloud1prime,5,'y',5)
%plotCloud(matchcloud2prime,5,'c',5)
axis manual
axis(scale)
title('possible matches')
subplot(4,1,4)
plotCloud(matchcloud1p,5,'r',5)
hold on
plotCloud(matchcloud2p,5,'b',5)
hold on
line([matchcloud1p(2,:);matchcloud2p(2,:)],[matchcloud1p(1,:);matchcloud2p(1,:)],-[matchcloud1p(3,:);matchcloud2p(3,:)])
axis manual
axis(scale)
title('possible matches windowed by normal')



figure(6)
plotCloud(matchcloud1,6,'r',5)
hold on
plotCloud(matchcloud2,6,'b',5)
plotCloud(matchcloud1p,6,'y',5)
plotCloud(matchcloud2p,6,'c',5)
%% now compare descriptors
figure(3); subplot(2,2,1); plot(abs(PDfeatures1.Descriptors(4,dist<ball)'-PDfeatures2.Descriptors(4,idxInsideRadius)')); hold on;
%plot(PDfeatures2.Descriptors(4,idxInsideRadius)','r');
title('feature nz')
subplot(2,2,3); plot(PDfeatures2.Descriptors(5,dist<ball)','b'); hold on; 
plot(PDfeatures2.Descriptors(5,idxInsideRadius)','r'); title('grazing angle'); 
subplot(2,2,2); plot(PDfeatures1.Descriptors(6,dist<ball)','b');hold on;plot(PDfeatures2.Descriptors(6,idxInsideRadius)','r');title('k1'); 
subplot(2,2,4); plot(PDfeatures1.Descriptors(7,dist<ball)'); hold on; plot(PDfeatures2.Descriptors(7,idxInsideRadius)','r');title('k2'); 

% scale axes
scale = [min(cloud'); max(cloud')];
scale = reshape(scale,1,[]);
% 
% matchDescriptors1 = PDfeatures1.Descriptors(:,dist<ball);
% matchDescriptors2 = PDfeatures2.Descriptors(:,idxInsideRadius);
% 
% surf(matchDescriptors1'*matchDescriptors2)
