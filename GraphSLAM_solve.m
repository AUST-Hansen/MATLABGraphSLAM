function [mu, Sigma] = GraphSLAM_solve(OmegaReduced,zetaReduced,Omega,zeta,c_i_t)

    Sigma = inv(OmegaReduced);
    %Sigma = 0;
    muPoses = OmegaReduced\zetaReduced;
    size(muPoses)
    size(zetaReduced)



    %mu = [muPoses; muMap];
    nPoses = length(muPoses)/6;
    nMapFeatures = (size(Omega,2) - size(OmegaReduced,2))/3;
    muMap = zeros(3*nMapFeatures,1);
    for jj = 1:nMapFeatures
        
        OmegaJJ = Omega(6*nPoses+jj*3-2:6*nPoses+jj*3,6*nPoses+jj*3-2:6*nPoses+jj*3);
        zetaJ = zeta(6*nPoses+jj*3-2:6*nPoses+jj*3);
        OmegaJTau = Omega(6*nPoses+jj*3-2:6*nPoses+jj*3,1:6*nPoses);
        
        muMap(jj*3-2:jj*3) = OmegaJJ\(zetaJ - OmegaJTau*muPoses(1:6*nPoses));
        
        
    end

    mu = [muPoses;muMap]
    fprintf('Solve\n')
