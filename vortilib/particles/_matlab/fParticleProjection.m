function [MeshValues v_p]=fParticleProjection(Part, Mesh, Algo)

    %% Project Particles into a grid

    if Algo.nDim==2
        v_p =zeros(Part.nPart,2); % 1: intensity, 2: volume
        v_p(:,1)=Part.Intensity; % Gamma
        v_p(:,2)=Part.Volume;
    elseif Algo.nDim==3
        v_p =zeros(Part.nPart,4); % 1:3: intensity, 4: volume
        v_p(:,1:3)=Part.Intensity; % Gamma
        v_p(:,4)=Part.Volume;
    end

%     MeshValues= interp_p2m(Part, Mesh, v_p, Algo);
%       MeshValues= interp_p2m_mex(Part.P,Part.nPart,Mesh.n,v_p,Algo.nDim,Algo.InterpolationKernel,Mesh.v1,Mesh.v2,Mesh.v3,Mesh.bRegular);
      MeshValues= interp_p2m(Part.P,Part.nPart,Mesh.n,v_p,Algo.nDim,Algo.InterpolationKernel,Mesh.v1,Mesh.v2,Mesh.v3,Mesh.bRegular);
