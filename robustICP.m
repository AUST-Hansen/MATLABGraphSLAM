function [TR, TT, ER] = robustICP(pointCloud1,pointCloud2,varargin)


    worstRejection = 0.2;
    twoDee = true;
    
    lenP = length(pointCloud1)
    lenQ = length(pointCloud2)
    
    smoothingFactor = 7;
    smooth1 = [smooth(pointCloud1(1,:),smoothingFactor)';smooth(pointCloud1(2,:),smoothingFactor)';smooth(pointCloud1(3,:),smoothingFactor)'];
    smooth2 = [smooth(pointCloud2(1,:),smoothingFactor)';smooth(pointCloud2(2,:),smoothingFactor)';smooth(pointCloud2(3,:),smoothingFactor)'];
    
    smooth1Prime = diff(smooth1,1,2);
    smooth2Prime = diff(smooth2,1,2);
    smooth1DPrime = diff(smooth1,2,2);
    smooth2DPrime = diff(smooth2,2,2);
    % curvature approximation
    k1 = abs(smooth1Prime(1,2:end).*smooth1DPrime(2,:) - smooth1Prime(2,2:end).*smooth1DPrime(1,:))./((smooth1Prime(1,2:end).^2 + smooth1Prime(2,2:end).^2).^1.5);
    k2 = abs(smooth2Prime(1,2:end).*smooth2DPrime(2,:) - smooth2Prime(2,2:end).*smooth2DPrime(1,:))./((smooth2Prime(1,2:end).^2 + smooth2Prime(2,2:end).^2).^1.5);
    % unweigh edgeweights
    w1 = 15*[0,k1,0];% w1(1:15) = 0; w1(end-15:end) = 0
    w2 = 15*[0,k2,0];% w2(1:15) = 0; w2(end-15:end) = 0
    w1(isnan(w1)) = 0;
    w2(isnan(w2))=0;
% Do ICP
            %[TR,TT,ER] = icp(pointCloud1,pointCloud2,20,'WorstRejection',worstRejection,'twoDee',twoDee,'Weight',@weighting);
            [TR,TT,ER] = icp(smooth1,smooth2,20,'WorstRejection',worstRejection,'twoDee',twoDee,'Weight',@weighting);
            ER
            
            
            
            
    function weights = weighting(indices)
        indices
        size(w1)
        size(w2)
        w1
        w2
        weights = w1(indices).*w2;
        
    end   
end