% cleanMeasurementsTestHarness.m

addpath ..

useDVL = false;
useImagenex = true;
useMultibeam = false;
skip = 3;
start = 1;
stop = 1000; %size(rangeData,2);

[numMeasurements, rangeMeasurements, c_i_t] = cleanMeasurements(sensor,rangeData(:,start:skip:stop),useMultibeam,imagenex,imagenexData(:,start:skip:stop),useImagenex,dvl,dvlData(start:skip:stop),useDVL);