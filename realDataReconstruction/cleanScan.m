function cleaned_scan = cleanScan(scan,minClusterSize)

   % takes scan in vehicle frame( x should all be zero)
   
   threshold = 2; % meters difference at 40 m range
   
   ranges = sqrt(scan(1,:).^2 + scan(2,:).^2);
   angles = atan2(-scan(3,:),-scan(2,:));
   
   dr = diff(ranges);
   regionlabels = zeros(size(ranges));
   
   region = 1;
   
   regionlabels(1) = 1;
   regioncounter = regionlabels;
   for ii = 2:size(ranges,2)
      if abs(dr(ii-1)) > threshold
          region = region+1;
      end
      regionlabels(ii) = region;
      regioncounter(region) = regioncounter(region) + 1;
       
   end
   
   indices = find(regioncounter>minClusterSize);
   
   goodindices = false(size(regionlabels));
   
   for ii = 1:size(indices,2)
       goodindices = goodindices | regionlabels == indices(ii);
   end
   
   cleaned_scan = scan(:,goodindices);
   
   

