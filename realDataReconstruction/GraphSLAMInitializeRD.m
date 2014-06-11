function [SLAMdata, otherData] = GraphSLAMInitializeRD(PD)

    xy = PD.Pos(1:2,:);
    uv = smooth(PD.DVL(1:2,:),20);
    psi = PD.Psi';
    bias = 0*psi;
    
    % Initialize slam quantities
    SLAMdata.mu = reshape([xy;psi;uv;bias],[],1);
    SLAMdata.c_i_t = 
    % Other quantities
    