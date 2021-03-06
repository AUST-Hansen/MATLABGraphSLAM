function initialMapEstimate = GraphSLAM_initializeMap(initialStateEstimate,rangeMeasurements,correspondences)

   N = size(initialStateEstimate,2);
   M = size(correspondences.c_i_t,1);
   initialMapEstimate = 0*rangeMeasurements;
   jMapFeature = 1;
   for ii = 1:N
       
       B_R_V = Euler2RotMat(0,0,initialStateEstimate(5,ii)); 
       x_Veh_bergframe = [initialStateEstimate(1:2,ii); 0];
       
       for jj = 1:M
          if(correspondences.c_i_t(jj,ii) <= 0) % no more valid returns at this timestep
              continue;
          else
              % extract measurement index
              measIdx = correspondences.c_i_t(jj,ii);
              % build measurement vector in body frame
              % measurement is [range bearing elevation]'
              meas = rangeMeasurements(:,measIdx);
              
             
              initialMapEstimate(:,jMapFeature) = x_Veh_bergframe + B_R_V*meas;
              jMapFeature = jMapFeature+1;
          end
       end
       
       
     
       
       
   end
  

