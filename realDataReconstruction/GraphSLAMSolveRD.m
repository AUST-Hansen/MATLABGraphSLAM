function SLAMdataOut = GraphSLAMSolveRD(SLAMdata,OtherData);

    Sigma = inv(SLAMdata.OmegaReduced);
    %Sigma = 0;
    muPoses = SLAMdata.OmegaReduced\SLAMdata.zetaReduced;

    fprintf('Solve\n')

    %mu = [muPoses; muMap];
    nPoses = SLAMdata.Nstates;
    nMapFeatures = (size(SLAMdata.Omega,2) - size(SLAMdata.OmegaReduced,2))/SLAMdata.mapDimension;
    muMap = zeros(SLAMdata.mapDimension*nMapFeatures,1);
    for jj = 1:nMapFeatures
        
        OmegaJJ = SLAMdata.Omega(SLAMdata.stateSize*nPoses+(jj-1)*SLAMdata.mapDimension+1:...
                        SLAMdata.stateSize*nPoses+jj*SLAMdata.mapDimension,...
                        SLAMdata.stateSize*nPoses+(jj-1)*SLAMdata.mapDimension+1:...
                        SLAMdata.stateSize*nPoses+jj*SLAMdata.mapDimension);
        zetaJ = SLAMdata.zeta(SLAMdata.stateSize*nPoses+(jj-1)*SLAMdata.mapDimension+1:...
                        SLAMdata.stateSize*nPoses+jj*SLAMdata.mapDimension);
        OmegaJTau = SLAMdata.Omega(SLAMdata.stateSize*nPoses+(jj-1)*SLAMdata.mapDimension+1:...
                          SLAMdata.stateSize*nPoses+jj*SLAMdata.mapDimension,...
                          1:SLAMdata.stateSize*nPoses);
        

        muMap((jj-1)*SLAMdata.mapDimension+1:jj*SLAMdata.mapDimension) = OmegaJJ\(zetaJ - OmegaJTau*muPoses(1:SLAMdata.stateSize*nPoses));
        
        if(mod(jj,500)==0)
        fprintf('%d of %d features processed\n',jj,nMapFeatures);
    end
        
    end

    mu = [muPoses;muMap]

    SLAMdataOut = SLAMdata;
    SLAMdataOut.mu = mu;
    SLAMdataOut.Sigma = Sigma;