function [OmegaTilde,zetaTilde] = GraphSLAM_reduce(timeStamps,stateSize,Omega,zeta)

   lastStateIdx = length(timeStamps)*stateSize;
   nMapFeatures = size(zeta(lastStateIdx+1:end),1)/3;
   
   OmegaTilde = Omega(1:lastStateIdx,1:lastStateIdx);
   zetaTilde = zeta(1:lastStateIdx);

   for jj = 1:nMapFeatures

       OmegaTilde = OmegaTilde - sparse(Omega(1:lastStateIdx,lastStateIdx+3*(jj-1) + 1:lastStateIdx+3*jj)*...
                                 (Omega(lastStateIdx+3*(jj-1) + 1:lastStateIdx+3*jj,lastStateIdx+3*(jj-1) + 1:lastStateIdx+3*jj)\...
                                 Omega(lastStateIdx+3*(jj-1) + 1:lastStateIdx+3*jj,1:lastStateIdx)));
       zetaTilde = zetaTilde - Omega(1:lastStateIdx,lastStateIdx+3*(jj-1) + 1:lastStateIdx+3*jj)*...
                                 (Omega(lastStateIdx+3*(jj-1) + 1:lastStateIdx+3*jj,lastStateIdx+3*(jj-1) + 1:lastStateIdx+3*jj)\...
                                 zeta(lastStateIdx+3*(jj-1) + 1:lastStateIdx+3*jj));
       Omega(lastStateIdx+3*(jj-1) + 1:lastStateIdx+3*jj,lastStateIdx+3*(jj-1) + 1:lastStateIdx+3*jj);   
       
       % Visualize progress
       if(mod(jj,100) == 0)
           fprintf('%d of %d\n',jj,nMapFeatures)
       end
   end
  