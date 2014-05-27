function DataOut = extractData(file)

    M = csvread(file,1,0);
    % Extract time    
    [timeStamps, IA, IC] = unique(M(:,1));
    timeSteps = timeStamps - timeStamps(1);
    % Extract position, heading, speed
    LatLonDepth = M(IA,2:4);
    Psi = M(IA,5)*pi/180.;
    Speed = M(IA,6)*.3048 ; % I think this is in feet/sec, so convert to m/s   
    [eastings,northings,zone] = deg2utm(LatLonDepth(:,1),LatLonDepth(:,2));
    XYZ = [northings, eastings, LatLonDepth(:,3)]' - repmat([northings(1);eastings(1); 0],1,size(northings,1));
    % build correspondence table
    nBeams = 512;
    c_i_t = zeros(nBeams,length(timeSteps));
    rangecounter = 1;
    for ii = 1:length(timeSteps)    
        range_t = M(IC==ii,7:9)';
        m = length(range_t);
        if (m>nBeams) % duplicate record
            c_i_t(1:m/2,ii) = (rangecounter:rangecounter+(m/2)-1)';
        else
            c_i_t(1:m,ii) = (rangecounter:rangecounter+m-1)';
        end
        rangecounter = rangecounter+m;
        if (mod(ii,100) == 0)
            fprintf('processed %3.1f of %3.1f scans\n',ii,length(timeSteps))
        end
    end
    % Extract range measurements
    rangeMeasurements = M(:,7:9)';
    % Pack up all data
    DataOut.timeSteps = timeSteps;
    DataOut.Pos = XYZ;
    DataOut.Psi = Psi;
    DataOut.Phi = 0*Psi;
    DataOut.Theta = 0*Psi;
    DataOut.Speed = Speed;
    DataOut.c_i_t = sparse(c_i_t);
    DataOut.rangeMeasurements = rangeMeasurements;
