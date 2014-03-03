% robustscanmatch test harness
close all

step = poseSkip;
for ii = 1:4:193
    close all;
idx1 = ii;
idx2 = ii+step;

pointCloud1 = rangeMeasurements(:,meas_ind(meas_ind(:,idx1)~=-17,idx1));
pointCloud2 = rangeMeasurements(:,meas_ind(meas_ind(:,idx2)~=-17,idx2));

[ii]

[TR, TT, ER, knnOut] = robustScanMatch(pointCloud1,pointCloud2,'DEBUG');

yaw = asin(-TR(1,2))
180/pi*yaw

Input = input('Press enter to continue\n','s');

end