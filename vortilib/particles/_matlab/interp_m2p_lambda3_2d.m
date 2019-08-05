function part_values= interp_m2p_lambda3_2d(Part,nPart,mesh,nval,v1,v2,n1,n2,bRegular)
% xBox: Origin of the grid
% dBox: Vector containing the length of a grid cell in the three direction
% part       (npart, 1:ndim)
% part_values(npart_values, 1:nval)
% mesh(1:nval,nx,ny)

if nargin==0
    nPart=4;
    nDim=2;
    nVal=3;
    Part=zeros(nPart,nDim);
    Part(1,:) = [0,0];
    Part(2,:) = [1.5,0.5];
    Part(3,:) = [2.0,1.0];
    Part(4,:) = [4,4];
    v1 = [0 ,1,2,3 ,4];
    v2 = [-2,0,2,4];
    n1=length(v1);
    n2=length(v2)
    [X,Y]=meshgrid(v1,v2);
    mesh=zeros(nVal,n1,n2);
    mesh(1,:,:)=1+X'*0;
    mesh(2,:,:)=2+Y'*0;
    mesh(3,:,:)=X'.*Y';
    bRegular=true;
    part_values= interp_m2p_lambda3_2d(Part,nPart,mesh,nVal,v1,v2,n1,n2,bRegular)
    pythondisp(squeeze(part_values),'v_p')
    return
end

Part=Part'; % We transpose because code below is fortran oriented

if(size(mesh,1)~=nval)
    disp('Mesh has wrong size')
%     kbd
end
if(size(Part,2)~=nPart)
    disp('Part has wrong size')
%     kbd
end
% M4' kernel needs four points
%         xp
% i1   i2   i3   i4
%  |    | .  |    |

part_values=zeros(nval,nPart);
ic=[0,0];
dc=[0,0];
xBox=[0,0];
dBox=[0,0];
xBox(1)=v1(1);
xBox(2)=v2(1);
dBox(1)=v1(2)-v1(1);
dBox(2)=v2(2)-v2(1);
for i=1:nPart
    % Coordinates and distance in grid index space (to nearest left grid point)
    if bRegular==1
        [ic(1), dc(1)]=fCoordRegularGrid(Part(1,i),xBox(1),dBox(1));
        [ic(2), dc(2)]=fCoordRegularGrid(Part(2,i),xBox(2),dBox(2));
    else
        [ic(1), dc(1)]=fCoordRectilinearGrid(Part(1,i),v1);
        [ic(2), dc(2)]=fCoordRectilinearGrid(Part(2,i),v2);
    end

    % Getting the M'4 kernel coefficients
    [ax1,ax2,ax3,ax4,i1,i2,i3,i4]=interp_coeff_lambda3(ic(1),dc(1),n1);
    [ay1,ay2,ay3,ay4,j1,j2,j3,j4]=interp_coeff_lambda3(ic(2),dc(2),n2);

    % !---------------------------------------------------------------------------
    % ! Combine 1D kernels into the 3D kernel and apply to the field
    % !---------------------------------------------------------------------------
    part_values(1:nval,i) = ax1*ay1*mesh(1:nval,i1,j1) + ...
                            ax2*ay1*mesh(1:nval,i2,j1) + ...
                            ax3*ay1*mesh(1:nval,i3,j1) + ...
                            ax4*ay1*mesh(1:nval,i4,j1) + ...
                            ax1*ay2*mesh(1:nval,i1,j2) + ...
                            ax2*ay2*mesh(1:nval,i2,j2) + ...
                            ax3*ay2*mesh(1:nval,i3,j2) + ...
                            ax4*ay2*mesh(1:nval,i4,j2) + ...
                            ax1*ay3*mesh(1:nval,i1,j3) + ...
                            ax2*ay3*mesh(1:nval,i2,j3) + ...
                            ax3*ay3*mesh(1:nval,i3,j3) + ...
                            ax4*ay3*mesh(1:nval,i4,j3) + ...
                            ax1*ay4*mesh(1:nval,i1,j4) + ...
                            ax2*ay4*mesh(1:nval,i2,j4) + ...
                            ax3*ay4*mesh(1:nval,i3,j4) + ...
                            ax4*ay4*mesh(1:nval,i4,j4) ;
end


part_values=part_values'; % We transpose because code above is fortran oriented
