function [OmegaTilde,zetaTilde] = GraphSLAM_reduce(timeStamps,stateSize,Omega,zeta,c_i_t)

   lastStateIdx = length(timeStamps)*stateSize;
   nMapFeatures = max(max(c_i_t)) ; %size(zeta(lastStateIdx+1:end),1)/3;
   
   OmegaTilde = Omega(1:lastStateIdx,1:lastStateIdx);
   zetaTilde = zeta(1:lastStateIdx);


   for jj = 1:nMapFeatures
       jj
       % identify poses which have seen map features
       [aa,bb]=find(c_i_t == jj);
       if (~isempty(bb)) % If we've seen this feature
       
           OmegaJJ = Omega(lastStateIdx+3*jj-2:lastStateIdx+3*jj,lastStateIdx+3*jj-2:lastStateIdx+3*jj); 
           OmegaTauJ = Omega(1:stateSize*size(c_i_t,2),lastStateIdx+3*jj-2:lastStateIdx+3*jj);
           OmegaTilde = OmegaTilde - sparse(OmegaTauJ*(OmegaJJ\(OmegaTauJ')));
           zetaTilde = zetaTilde - OmegaTauJ*(OmegaJJ\zeta(lastStateIdx+3*jj-2:lastStateIdx+3*jj));
       end
   end
   % Visualize progress
   if(mod(jj,100) == 0)
       fprintf('%d of %d\n',jj,nMapFeatures)
   end
   
end
