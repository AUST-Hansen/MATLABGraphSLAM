function mu_out = GraphSLAM_initialize(timeStamps,inputHistory,initialPsiDotGuess)

   % inputHistory(t) = [vx_nominal omega_in_worldframe]'
   N = length(inputHistory);
   
   % State: [x_bergframe y_bergframe u v heading_bergframe omega_berg]'
   mu_out = zeros(6,N);
   mu_out(3,1) = inputHistory(1,1);
   mu_out(6,1) = initialPsiDotGuess;
   
   for ii = 2:N
       dT = timeStamps(ii) - timeStamps(ii-1);
       heading = mu_out(5,ii-1);
       berg_R_veh = [cos(heading), -sin(heading); sin(heading), cos(heading)];
       mu_out(1:2,ii) = eye(2)*mu_out(1:2,ii-1) + berg_R_veh*mu_out(3:4,ii-1)*dT; % update position
       mu_out(3:4,ii) = .2*mu_out(3:4,ii-1) + .8*[inputHistory(1,ii); 0];
       mu_out(5,ii) = mu_out(5,ii-1) + (inputHistory(2,ii) - mu_out(6,ii-1))*dT;
       mu_out(6,ii) = mu_out(6,ii-1); 
   end
   
   % done