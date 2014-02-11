% robustscanmatch test harness
close all

idx1 = 9;
idx2 = 13;
pointCloud1 = rangeMeasurements(:,c_i_t(c_i_t(:,idx1)~=-17,idx1));
pointCloud2 = rangeMeasurements(:,c_i_t(c_i_t(:,idx2)~=-17,idx2));



[TR, TT, ER] = robustScanMatch(pointCloud1,pointCloud2);