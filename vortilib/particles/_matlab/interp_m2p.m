function v_p= interp_m2p(PartP,nPart,n,MeshValues,nDim,InterpolationKernel,v1,v2,v3,bRegular);

    % Some simple variables
%     nDim=Algo.nDim;
%     xBox=Mesh.xMesh_min;
%     dBox=Mesh.dCell;
%     n=Mesh.n;
    
%     P=Part.P; %  nPart x nDim
% %     P=checkNcolumns(P,nDim,'P');
%     nPart=Part.nPart;
    nVal=size(MeshValues,1); % MeshValues is nVal * nx * ny for Now

    v_p=0;
    % Since matlab coder dont like switch with strings
    if isequal(lower(InterpolationKernel),'mp4')
        s=1;
    else
        s=2; %lambda 3
    end

    switch(s)
        case 1
           switch(nDim) 
               case 2
                   v_p= interp_m2p_mp4_2d(PartP,nPart,MeshValues,nVal,v1,v2,n(1),n(2),bRegular);
%                    check_size(v_p,[nPart nVal],'v_p'); 
           end
        case 2
           switch(nDim) 
               case 2
                   v_p= interp_m2p_lambda3_2d(PartP,nPart,MeshValues,nVal,v1,v2,n(1),n(2),bRegular);
%                    check_size(v_p,[nPart nVal],'v_p'); 
           end


        otherwise
%             log_error('Unknown Interpolation kernel')

    end


