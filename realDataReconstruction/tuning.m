% tuning script
close all;
%%
PDfeatures = sweepSonarFeatures(ProcessedData);
%%
cloud = buildFeatureMapFromProcessedData(PDfeatures,1,size(PDfeatures.c_i_t,2));
cloud2 = buildMapFromProcessedData(ProcessedData,1,size(ProcessedData.c_i_t,2));
plotCloud(cloud2(:,40000:100:end),5,'r');
hold on;
plotCloud(cloud,5,'g')
%%
figure(3); plot(PDfeatures.Descriptors(4:end,:)'),legend('nz','grazing angle','k1','k2')