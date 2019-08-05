function [Part, Field]=fParticleInitialization(Part, InitParams, Algo)

    Field=[];
    switch(InitParams.Method)

        case('Mesh_GridPoints_OmegaAnalytical')
            Mesh=InitParams.Mesh;

            % Putting particles at the middle of the cell
            cell_volumes = Mesh.getCellVolumes()     ;  % This is cell area in 2D !!!!
            PP           = Mesh.getFlatGridPoints() ; 
            nPart=size(PP,1);
            Part.reset(nPart);
            Part.setP(PP);

            % Computing Omega from analytical function and setting particle intensities
            if Algo.nDim==2
                Omega=InitParams.fOmega(Part.P(:,1),Part.P(:,2));
            else
                Omega=InitParams.fOmega(Part.P(:,1),Part.P(:,2),Part.P(:,3));
            end
            Part.setIntensity(Omega.*cell_volumes); % alpha_p=omega*V  ; Gamma=omega*A
            Part.setVolume(cell_volumes);

        case('Mesh_CellCenter_OmegaAnalytical')
            Mesh=InitParams.Mesh;

            % Putting particles at the middle of the cell
            cell_volumes=Mesh.getCellVolumes();  % This is cell area in 2D !!!!
            cell_Centers=Mesh.getFlatCellCenters();
            nPart=size(cell_Centers,1);
            Part.reset(nPart);
            Part.setP(cell_Centers);

            % Computing Omega from analytical function and setting particle intensities
            if Algo.nDim==2
                Omega=InitParams.fOmega(Part.P(:,1),Part.P(:,2));
            else
                Omega=InitParams.fOmega(Part.P(:,1),Part.P(:,2),Part.P(:,3));
            end
            Part.setIntensity(Omega.*cell_volumes); % alpha_p=omega*V  ; Gamma=omega*A
            Part.setVolume(cell_volumes);
        
        case('Mesh_CellCenterQuasiRandom_OmegaAnalytical')
            Mesh=InitParams.Mesh;
            if ~Mesh.bRegular
                error('Random init with non regular grid not suported')
            end

            
            % Putting particles at the middle of the cell
            cell_volume=Mesh.getCellVolume();  % This is cell area in 2D !!!!
            cell_Centers=Mesh.getFlatCellCenters();
            nPart=size(cell_Centers,1);

            % rand is between [0 1] we transform it between -1/2h 1/2h
            RandP=zeros(nPart,Algo.nDim);
            for id=1:Algo.nDim
                RandP(:,id)=(rand(nPart,1)-1/2)*Mesh.dCell(id);
            end

            % Setting Particles position
            Part.reset(nPart);
            Part.setP(cell_Centers+RandP);

            % Computing Omega from analytical function and setting particle intensities
            if Algo.nDim==2
                Omega=InitParams.fOmega(Part.P(:,1),Part.P(:,2));
            else
                Omega=InitParams.fOmega(Part.P(:,1),Part.P(:,2),Part.P(:,3));
            end
            Part.setIntensity(Omega*cell_volume); % alpha_p=omega*V  ; Gamma=omega*A
            Part.setVolume(cell_volume);
        
        case('Mesh_Random_OmegaAnalytical')
            Mesh=InitParams.Mesh;

            
            % finding spatial extent and total number
            nPart=1;
            TotalVolume=1;
            extent=zeros(1,Algo.nDim);
            xmin=zeros(1,Algo.nDim);
            for id=1:Algo.nDim
                nPart=nPart*Mesh.n(id);
                extent(id)=Mesh.xMesh_max(id)-Mesh.xMesh_min(id);
                xmin(id)=Mesh.xMesh_min(id);
                TotalVolume=TotalVolume*extent(id);
            end

            part_volume=TotalVolume/nPart; 


            % rand is between [0 1] we transform it between  [0 extent] + xmin
            RandP=zeros(nPart,Algo.nDim);
            for id=1:Algo.nDim
                RandP(:,id)=rand(nPart,1)*extent(id) + xmin(id);
            end

            %  Setting Particle Positions
            Part.reset(nPart);
            Part.setP(RandP);

            % Computing Omega from analytical function and setting particle intensities
            if Algo.nDim==2
                Omega=InitParams.fOmega(Part.P(:,1),Part.P(:,2));
            else
                Omega=InitParams.fOmega(Part.P(:,1),Part.P(:,2),Part.P(:,3));
            end
            Part.setIntensity(Omega*part_volume); % alpha_p=omega*V  ; Gamma=omega*A
            Part.setVolume(part_volume);

        case('Mesh_CellCenter_OmegaM2P')
            Mesh=InitParams.Mesh;

            % Putting particles at the middle of the cell
            cell_volumes =Mesh.getCellVolumes();  % This is cell area in 2D !!!!
            cell_Centers =Mesh.getFlatCellCenters();
            nPart=size(cell_Centers,1);
            Part.reset(nPart);
            Part.setP(cell_Centers);

            v_p= interp_m2p_mex(Part.P,Part.nPart,Mesh.n, InitParams.MeshValues, Algo.nDim, Algo.InterpolationKernel,Mesh.v1,Mesh.v2,Mesh.v3,Mesh.bRegular);
            if Algo.nDim==2
                Part.setIntensity(v_p(:,1)); %  Gamma=omega*A
            else
                Part.setIntensity(v_p(:,1:3)); % alpha_p=omega*V 
            end
            Part.setVolume(cell_volumes); 

        case('Mesh_GridPoints_OmegaM')
            Mesh=InitParams.Mesh;

            % Putting particles at the middle of the cell
            cell_volumes=Mesh.getFakeCellVolumes();  % This is cell area in 2D !!!!
            grid_points =Mesh.getFlatGridPoints();
            nPart=size(grid_points,1);
            Part.reset(nPart);
            Part.setP(grid_points);

%             v_p= interp_m2p_mex(Part.P,Part.nPart,Mesh.xMesh_min,Mesh.dCell,Mesh.n, InitParams.MeshValues, Algo.nDim, Algo.InterpolationKernel);
            MeshVal_flat=Mesh.flattenValues(InitParams.MeshValues);
            if Algo.nDim==2
                MeshIntensity = MeshVal_flat(1,:)' ; 
                MeshVol       = MeshVal_flat(2,:)' ; 
            else
                MeshIntensity = MeshVal_flat(1:3,:)' ; 
                MeshVol       = MeshVal_flat(2,:)' ; 
            end
%             kbd
%             MeshVol(MeshVol<=0)=cell_volumes(MeshVol<=0);
            MeshOmega=MeshIntensity./MeshVol;
%             kbd
            Part.setIntensity(MeshOmega(:).*cell_volumes(:)); % alpha_p=omega*V 
%             Part.setIntensity(MeshIntensity); % alpha_p=omega*V 
            Part.setVolume(cell_volumes); 

%             v_p= interp_m2p_mex(Part.P,Part.nPart,Mesh.xMesh_min,Mesh.dCell,Mesh.n, InitParams.MeshValues, Algo.nDim, Algo.InterpolationKernel);
%             Part.setIntensity(v_p(:,1)); %  Gamma=omega*A



        otherwise
            log_error('Not implemented')
    end


end
