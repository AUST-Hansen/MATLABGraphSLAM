function [mu, Sigma] = GraphSLAM_solve(OmegaReduced,zetaReduced,Omega,zeta)


    Sigma = inv(OmegaReduced);
    
    % MATLAB seems to handle this pseudoinversion just fine
    % If it has problems later, do the procedure on p 349 of ProbRob
    mu = Omega\zeta;