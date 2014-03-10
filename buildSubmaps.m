function Submaps = buildSubmaps(Velocities,Eulers,Length,measurements, correspondences)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Function buildSubmaps
%
%   Inputs:
%      Velocities - 3xN matrix of body-frame velocities
%      Eulers - 3xN matrix of vehicle Euler angles
%      Length - scalar, number of timesteps to be concatenated
%      measurements - multibeam measurements
%
%   Outputs:
%      Submaps - array of submap structures containing:
%         Reference pose index (submaps are associated with the vehicle
%         pose at the 'middle' of the submap)
%         Point cloud of map points in the reference pose frame
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Submaps = [];
NumPings = size(Velocities,2);
NumSubmaps = floor(NumPings/Length);


%Build dead-reckoned trajectory
x_veh = 0*Velocities;
for ii = 2:NumPings
    x_veh(:,ii) = x_veh(:,ii-1) + Euler2RotMat(Eulers(1,ii),Eulers(2,ii),Eulers(3,ii))*Velocities(:,ii-1);
end

t = 1;
refInd = ceil(Length/2);
Submaps = [];
Submap = [];
for i_map = 1:NumSubmaps
    Subtrajectory = x_veh(:,i_map*Length-(Length-1):i_map*Length);
    SubEuler = Eulers(:,i_map*Length-(Length-1):i_map*Length);
    SubCIT = correspondences.c_i_t(:,i_map*Length-(Length-1):i_map*Length);
    Submap.refIndex = (i_map-1)*Length + refInd;
    I_R_ref = Euler2RotMat(SubEuler(1,refInd),SubEuler(2,refInd),SubEuler(3,refInd));
    cloud = [];
    % build up submap in inertial frame
    for i_ping = 1:Length
       submeas_ind = SubCIT(SubCIT(:,i_ping)>0,i_ping);
       I_R_V = Euler2RotMat(SubEuler(1,i_ping),SubEuler(2,i_ping),SubEuler(3,i_ping));
       if (~isempty(submeas_ind))
           cloud = [cloud, (I_R_V*measurements(:,submeas_ind) + repmat(Subtrajectory(:,i_ping),1,length(submeas_ind)))];
       end
    end
    % now put into reference frame of reference index
    Submap.cloud = I_R_ref'*(cloud - repmat(Subtrajectory(:,refInd),1,size(cloud,2)));
    % push onto the stack
    Submaps = [Submaps, Submap];
end