function [pointCloud, x_Veh_bergframe] = reconstructTrueCloud(vehiclePosition,vehicleAttitude,icebergPosition,icebergAttitude,sensor,rangeData,varargin)
%--------------------------------------------------------------------------
% Reconstruct either imagenex or reson data from truth
%  Inputs:
%       vehiclePosition = 3xN matrix of xyz in inertial 
%       vehicleHeading = 
%
%

   N = size(vehiclePosition,2);
   pointCloud = zeros(3,size(rangeData,2));
   iValidPoint = 0;
   x_Veh_bergframe = zeros(3,N);
   for ii = 1:N % for each pose
    
    % Rotation matrix from vehicle to berg
    I_R_B = Euler2RotMat(icebergAttitude(1,ii),icebergAttitude(2,ii),icebergAttitude(3,ii));
    I_R_V = Euler2RotMat(vehicleAttitude(1,ii),vehicleAttitude(2,ii),vehicleAttitude(3,ii));
    B_R_V = I_R_B'*I_R_V; 
    x_Veh_bergframe(:,ii) = I_R_B'*(vehiclePosition(:,ii) - icebergPosition(:,ii));
    
    for jj = 1:sensor.numBeams
       % rotation from beam frame to sensor frame
       S_R_B       = [ cos(sensor.az(jj))  , -sin(sensor.az(jj)), 0   ;
                    sin(sensor.az(jj))  , cos(sensor.az(jj)), 0   ;
                    0                   , 0                 , 1  ];

       if (rangeData(jj,ii) ~= -17 && ~isnan(rangeData(jj,ii)) && rangeData(jj,ii) < sensor.maxRange )    
           pointReturn = x_Veh_bergframe(:,ii)  + B_R_V*sensor.Veh_R_sen*S_R_B*[rangeData(jj,ii) 0 0]';
           iValidPoint = iValidPoint+1;
           pointCloud(:,iValidPoint) = pointReturn;
       end    
    
    end
    
    
    
   end
   
   pointCloud = pointCloud(:,1:iValidPoint);
   fprintf('%i points!\n',iValidPoint)
   % handle options
   i_vararg = 1;
   while (i_vararg <= nargin - 6)
     if( strcmp(varargin{i_vararg},'plot'))
         fprintf('plotting data\n')
         sparsity = 50;
         figure;
         scatter3(pointCloud(1,1:sparsity:end),pointCloud(2,1:sparsity:end),pointCloud(3,1:sparsity:end),ones(1,size(pointCloud(1,1:sparsity:end),2)),'b')
         view(0,90)
         axis equal
         title('true reconstruction')
         i_vararg = i_vararg+1;
     end
   end