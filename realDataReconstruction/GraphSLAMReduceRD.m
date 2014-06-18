function SLAMdataOut = GraphSLAMReduceRD(SLAMdata,OtherData);


   lastStateIdx = SLAMdata.Nstates*SLAMdata.stateSize;
   nMapFeatures = max(max(SLAMdata.c_i_t)) ; %size(zeta(lastStateIdx+1:end),1)/3;
   
   OmegaTilde = SLAMdata.Omega(1:lastStateIdx,1:lastStateIdx);
   zetaTilde = SLAMdata.zeta(1:lastStateIdx);

fprintf('Reducing...\n');
for jj = 1:nMapFeatures
    if(mod(jj,500)==0)
        fprintf('%d of %d features processed\n',jj,nMapFeatures);
    end
    % identify poses which have seen map features
    [aa,bb]=find(SLAMdata.c_i_t == jj);
    if (~isempty(bb)) % If we've seen this feature
        try
            OmegaJJ = SLAMdata.Omega(lastStateIdx+SLAMdata.mapDimension*(jj-1)+1:lastStateIdx+SLAMdata.mapDimension*jj,...
                lastStateIdx+SLAMdata.mapDimension*(jj-1)+1:lastStateIdx+SLAMdata.mapDimension*jj);
            OmegaTauJ = SLAMdata.Omega(1:SLAMdata.stateSize*size(SLAMdata.c_i_t,2),...
                lastStateIdx+SLAMdata.mapDimension*(jj-1)+1:lastStateIdx+SLAMdata.mapDimension*jj);
            OmegaTilde = OmegaTilde - sparse(OmegaTauJ*(OmegaJJ\(OmegaTauJ')));
            zetaTilde = zetaTilde - ...
                OmegaTauJ*(OmegaJJ\SLAMdata.zeta(lastStateIdx+SLAMdata.mapDimension*(jj-1)+1:lastStateIdx+SLAMdata.mapDimension*jj));
        catch uhoh
            keyboard
        end
    end

end

   SLAMdataOut = SLAMdata; 
   SLAMdataOut.OmegaReduced = OmegaTilde;
   SLAMdataOut.zetaReduced = zetaTilde;
   
end