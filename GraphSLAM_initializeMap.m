function initialMapEstimate = GraphSLAM_initializeMap(initialStateEstimate,rangeMeasurements,c_i_t)

   N = size(initialStateEstimate,2);
   M = size(c_i_t,1);
   initialMapEstimate = 0*rangeMeasurements;
   jMapFeature = 1;
   for ii = 1:N
       
       B_R_V = Euler2RotMat(0,0,initialStateEstimate(5,ii)); 
       x_Veh_bergframe = [initialStateEstimate(1:2,ii); 0];
       
       for jj = 1:M
          if(c_i_t(jj,ii) <= 0) % no more valid returns at this timestep
              continue;
          else
              % extract measurement index
              measIdx = c_i_t(jj,ii);
              % build measurement vector in body frame
              % measurement is [range bearing elevation]'
              meas = rangeMeasurements(:,measIdx);
              measVec = Euler2RotMat(0,0,meas(2))*Euler2RotMat(0,meas(3),0)*[meas(1);0;0];
             
              initialMapEstimate(:,jMapFeature) = x_Veh_bergframe + B_R_V*measVec;
              jMapFeature = jMapFeature+1;
          end
       end
       
       
     
       
       
   end
  

