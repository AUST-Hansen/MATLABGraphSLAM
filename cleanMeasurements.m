function [measurementTimestamps, rangeMeasurements, correspondences,meas_ind,multibeamMeasurements,correspondences_mb] = cleanMeasurements(timeSteps,multibeam,multibeamData,useMultibeam,imagenex,imagenexData,useImagenex,dvl,dvlData,useDVL,rangeSkipReson,rangeSkipImagenex,poseSkip)

   % rangeSkip gives density in beams
   % poseSkip gives how long to wait between imagenex for matching.
   N = length(dvlData);
   % each column of rangeMeasurements is [range bearing elevation]' 
   rangeMeasurements = zeros(3,N*(multibeam.numBeams + imagenex.numBeams + dvl.numBeams));
   multibeamMeasurements = rangeMeasurements;
   correspondences.c_i_t = -17*ones((multibeam.numBeams + imagenex.numBeams + dvl.numBeams),N);
   correspondences.weights = 0*correspondences.c_i_t;
   correspondences_mb.c_i_t = correspondences.c_i_t;
   correspondences.weights = 0*correspondences.c_i_t;
   j_validMeas = 0;
   j_validMBMeas = 0;
   % for each timestep
   for i_time = 1:poseSkip:N
       correspondenceCounter = 1;
       mbCounter = 1;
       if(useImagenex)
       for jj = 1:rangeSkipImagenex:imagenex.numBeams
                  % rotation from beam frame to sensor frame
          S_R_B       = [ cos(imagenex.az(jj))  , -sin(imagenex.az(jj)), 0   ;
                          sin(imagenex.az(jj))  , cos(imagenex.az(jj)), 0   ;
                          0                   , 0                 , 1  ];
          if (imagenexData(jj,i_time) ~= -17 && ~isnan(imagenexData(jj,i_time)) && imagenexData(jj,i_time) < imagenex.maxRange)
              % valid reading!
              j_validMeas = j_validMeas+1;
              % rotation matrix
              v_R_beam = imagenex.Veh_R_sen*S_R_B;
              % rotate into vehicle frame
              beamInVehicleFrame = v_R_beam(:,1)*imagenexData(jj,i_time);
              %[r,b,elev] = getRBE(beamInVehicleFrame);
              % changing measurement model to straight vector
              %rangeMeasurements(:,j_validMeas) = [r,b,elev]';
              rangeMeasurements(:,j_validMeas) = beamInVehicleFrame;
              correspondences.c_i_t(correspondenceCounter,i_time) = j_validMeas;
              correspondenceCounter = correspondenceCounter+1;
          end
       end
       end
       
       
       % check multibeam
       if (useMultibeam)
           for jj = 1:rangeSkipReson:multibeam.numBeams
               % rotation from beam frame to sensor frame
               S_R_B       = [ cos(multibeam.az(jj))  , -sin(multibeam.az(jj)), 0   ;
                   sin(multibeam.az(jj))  , cos(multibeam.az(jj)), 0   ;
                   0                   , 0                 , 1  ];
               if (multibeamData(jj,i_time) ~= -17 && ~isnan(multibeamData(jj,i_time)) &&  multibeamData(jj,i_time) < multibeam.maxRange)
                   % valid reading!
                   j_validMeas = j_validMeas+1;
                   % rotation matrix
                   v_R_beam = multibeam.Veh_R_sen*S_R_B;
                   % rotate into vehicle frame
                   beamInVehicleFrame = v_R_beam(:,1)*multibeamData(jj,i_time);
                   [r,b,elev] = getRBE(beamInVehicleFrame);
                   %rangeMeasurements(:,j_validMeas) = [r,b,elev]';
                   rangeMeasurements(:,j_validMeas) = beamInVehicleFrame;
                   correspondences.c_i_t(correspondenceCounter,i_time) = j_validMeas;
                   correspondenceCounter = correspondenceCounter+1;
               end
           end
       end
       % loop closure multibeam-only stuff
       for jj = 1:rangeSkipReson:multibeam.numBeams
           % rotation from beam frame to sensor frame
           S_R_B       = [ cos(multibeam.az(jj))  , -sin(multibeam.az(jj)), 0   ;
               sin(multibeam.az(jj))  , cos(multibeam.az(jj)), 0   ;
               0                   , 0                 , 1  ];
           if (multibeamData(jj,i_time) ~= -17 && ~isnan(multibeamData(jj,i_time)) &&  multibeamData(jj,i_time) < multibeam.maxRange)
               % valid reading!
               j_validMBMeas = j_validMBMeas+1;
               % rotation matrix
               v_R_beam = multibeam.Veh_R_sen*S_R_B;
               % rotate into vehicle frame
               beamInVehicleFrame = v_R_beam(:,1)*multibeamData(jj,i_time);
               [r,b,elev] = getRBE(beamInVehicleFrame);
               %rangeMeasurements(:,j_validMeas) = [r,b,elev]';
               multibeamMeasurements(:,j_validMBMeas) = beamInVehicleFrame;
               correspondences_mb.c_i_t(mbCounter,i_time) = j_validMBMeas;
               mbCounter = mbCounter+1;
           end
       end
       
       if (useDVL)
           for jj = 1:dvl.numBeams
               % rotation from beam frame to sensor frame
               
               if (dvlData(i_time).ranges(jj) ~= -17 && ~isnan(dvlData(i_time).ranges(jj)) && dvlData(i_time).ranges(jj) < dvl.maxRange)
                   % valid reading!
                   j_validMeas = j_validMeas+1;
                   % rotate into vehicle frame
                   beamInVehicleFrame = (dvl.Veh_R_dvl*dvl.beamsVF(:,jj))*dvlData(i_time).ranges(jj);
                   %[r,b,elev] = getRBE(beamInVehicleFrame);
                   %rangeMeasurements(:,j_validMeas) = [r,b,elev]';
                   rangeMeasurements(:,j_validMeas) = beamInVehicleFrame;
                   correspondences.c_i_t(correspondenceCounter,i_time) = j_validMeas;
                   correspondenceCounter = correspondenceCounter+1;
               end
           end
       end
   end
   % trim rangeMeasurements of unfilled entries
   rangeMeasurements = rangeMeasurements(:,1:j_validMeas);
   multibeamMeasurements = multibeamMeasurements(:,1:j_validMBMeas);
   measurementTimestamps = timeSteps;
   meas_ind = correspondences.c_i_t;
end


function [range, bearing, elevation] = getRBE(vector)

range = norm(vector);
bearing = atan2(vector(2),vector(1));
if (bearing < 0) % want bearings to be compass bearings
    bearing = bearing+2*pi;
end
elevation = asin(-vector(3)/(range+eps));

end