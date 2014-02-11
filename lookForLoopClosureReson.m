function LCobjs = lookForLoopClosureReson(Submaps)

DEBUG = true;
ERRs = 100;
LCobjs = [];

for iRef = 1:length(Submaps)
        ERRs = 100;
    for ii=iRef+1:length(Submaps)
        [R, T,ERR]=icp(Submaps(iRef).cloud,Submaps(ii).cloud,20,'WorstRejection',.1);
        ERRs = [ERRs, ERR(end)];
        
    end
    
    [y, i] = min(ERRs);
    
    if (DEBUG)
        ii = i + (iRef - 1);
        [R, T,ERR]=icp(Submaps(iRef).cloud,Submaps(ii).cloud,20,'WorstRejection',.1);
        c2 = R*Submaps(ii).cloud + repmat(T,1,size(Submaps(ii).cloud,2));
        scatter3(Submaps(iRef).cloud(1,:),Submaps(iRef).cloud(2,:),Submaps(iRef).cloud(3,:),'b');
        hold on
        scatter3(c2(1,:),c2(2,:),c2(3,:),'r');
        axis equal
        
        iRef
        ii
        ERRs
        %keyboard
        hold off;
        
        Input = input('Does this look good? Return = yes, q = quit, anything else = no ','s');
    end
    if(strcmp(Input,'q'))
        return
    end
    if (isempty(Input))
        LCobj.idx1 = Submaps(iRef).refIndex;
        LCobj.idx2 = Submaps(i).refIndex;
        LCobj.R_1_2 = R;
        LCobj.T_1_2 = T;
        
        LCobjs = [LCobjs LCobj]
    end
end