function MeshValues= interp_p2m(Part,nPart,n,v_p,nDim,InterpolationKernel,v1,v2,v3,bRegular);
% function v_p= interp_m2p(PartP,nPart,n,MeshValues,nDim,InterpolationKernel,v1,v2,v3,bRegular);


    % Some simple variables
%     nDim=Algo.nDim;
%     xBox=Mesh.xMesh_min;
%     dBox=Mesh.dCell;
%     v1=Mesh.v1;
%     v2=Mesh.v2;
%     n=Mesh.n;
% 
%     P=Part.P; %  nPart x nDim
%     P=checkNcolumns(P,nDim,'P');
%     nPart=Part.nPart;
    nVal=size(v_p,2); % v_p is nPart x nVal
    MeshValues=0;

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
                    MeshValues=interp_p2m_mp4_2d(Part,nPart,v_p,nVal,v1,v2,n(1),n(2),bRegular);
%                     check_size(MeshValues,[nVal n(1) n(2)],'MeshValues'); % We are keeping fortran convention for now
            end
        case 2
            switch(nDim) 
                case 2
                    MeshValues=interp_p2m_lambda3_2d(Part,nPart,v_p,nVal,v1,v2,n(1),n(2),bRegular);
%                    check_size(MeshValues,[nVal n(1) n(2)],'MeshValues'); % We are keeping fortran convention for now
            end

        otherwise
%             log_error('Unknown Interpolation kernel')
            disp('Unknown Interpolation kernel')
            MeshValues=0;

    end

