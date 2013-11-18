function [mu, Sigma] = GraphSLAM_solve(OmegaReduced,zetaReduced,Omega,zeta,c_i_t)

    %Sigma = inv(OmegaReduced);
    Sigma = 0;
    %muPoses = OmegaReduced\zetaReduced;
    % MATLAB seems to handle this pseudoinversion just fine
    % If it has problems later, do the procedure on p 349 of ProbRob
%     numMapFeatures = (size(Omega,1) - size(OmegaReduced,1))/3;
%     nPoses = length(zetaReduced)/6;
%     
%     muMap = zeros(numMapFeatures*3,1);
%     for jj = 1:numMapFeatures
% 
%         omegaJ = Omega(nPoses*6+3*jj-2:nPoses*6+3*jj,nPoses*6+3*jj-2:nPoses*6+3*jj);
% 
%         muJ = omegaJ\zeta(nPoses*6+3*jj-2:nPoses*6+3*jj);
%         
%         [~,poseIdx] = find(c_i_t == jj);
%         for kk = 1:length(poseIdx)
%             omegaPJ = Omega(nPoses*6+3*jj-2:nPoses*6+3*jj, 6*kk-5:6*kk);
%             muJ = muJ + omegaJ\(omegaPJ*muPoses(6*kk-5:6*kk));
%         end
%         muMap(3*jj-2:3*jj) = muJ;
%     end
%   
    %mu = [muPoses; muMap];
    
    mu = Omega\zeta;