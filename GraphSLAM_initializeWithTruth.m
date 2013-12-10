function mu_out = GraphSLAM_initializeWithTruth(timeStamps,vehicleHeading,icebergHeading,inertialPos,bergPos)

   % inputHistory(t) = [vx_nominal omega_in_worldframe]'
   N = length(timeStamps);
   
   % State: [x_bergframe y_bergframe u v heading_bergframe omega_berg]'
   mu_out = zeros(6,N);
   mu_out(3,1) = 2;
   mu_out(6,1) = 0;
   
   for ii = 2:N
       
       I_R_B = Euler2RotMat(0,0,icebergHeading(3,ii));
       I_R_V = Euler2RotMat(0,0,vehicleHeading(3,ii));
       B_R_V = I_R_B'*I_R_V;
       
       dT = timeStamps(ii) - timeStamps(ii-1);

       mu_out(1:2,ii) = B_R_V(1:2,1:2)*(inertialPos(1:2,ii) - bergPos(1:2,ii)) - B_R_V(1:2,1:2)*(inertialPos(1:2,1) - bergPos(1:2,1));
       mu_out(3:4,ii) = B_R_V(1:2,1:2)'*((mu_out(1:2,ii) - mu_out(1:2,ii-1))/dT);
       mu_out(5,ii) = vehicleHeading(3,ii) - icebergHeading(3,ii);
       mu_out(6,ii) = mu_out(6,ii-1); 
   end
   

   mu_out(6,1:end-1) = diff(icebergHeading(3,1:N))/dT;
  
   
   % done