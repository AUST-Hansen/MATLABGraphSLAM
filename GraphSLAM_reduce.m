function [OmegaTilde,zetaTilde] = GraphSLAM_reduce(timeStamps,stateSize,Omega,zeta,correspondences)

   lastStateIdx = length(timeStamps)*stateSize;
   nMapFeatures = max(max(correspondences.c_i_t)) ; %size(zeta(lastStateIdx+1:end),1)/3;
   
   OmegaTilde = Omega(1:lastStateIdx,1:lastStateIdx);
   zetaTilde = zeta(1:lastStateIdx);

fprintf('Reducing...\n');
   for jj = 1:nMapFeatures
       if(mod(jj,1000)==0)
           fprintf('%d of %d features processed\n',jj,nMapFeatures);
       end
       % identify poses which have seen map features
       [aa,bb]=find(correspondences.c_i_t == jj);
       if (~isempty(bb)) % If we've seen this feature
           try
           OmegaJJ = Omega(lastStateIdx+3*jj-2:lastStateIdx+3*jj,lastStateIdx+3*jj-2:lastStateIdx+3*jj); 
           OmegaTauJ = Omega(1:stateSize*size(correspondences.c_i_t,2),lastStateIdx+3*jj-2:lastStateIdx+3*jj);
           OmegaTilde = OmegaTilde - sparse(OmegaTauJ*(OmegaJJ\(OmegaTauJ')));
           zetaTilde = zetaTilde - OmegaTauJ*(OmegaJJ\zeta(lastStateIdx+3*jj-2:lastStateIdx+3*jj));
           catch
               keyboard
           end
       end
   end
   % Visualize progress
   if(mod(jj,100) == 0)
       fprintf('%d of %d\n',jj,nMapFeatures)
   end
   
end
