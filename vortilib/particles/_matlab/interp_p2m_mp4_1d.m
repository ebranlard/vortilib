function mesh= p2m_mp4_1d(Part,part_values,xBox,dBox,nx)
% xBox: Origin of the grid
% dBox: Vector containing the length of a grid cell in the three direction
% part_values(1:3, npart_values)
% part_values(nd, npart_values)

nPart=size(Part,2);
nd=size(part_values,1);


if(nx<4)
    disp('No guarantee with tiny grid ')
end


% M4' kernel needs four points
%         xp
% i1   i2   i3   i4
%  |    | .  |    |

mesh=zeros(nd,nx);
for i=1:nPart
    C(1) = (Part(1,i)-xBox(1))/dBox(1); % position in grid coordinates
    % Getting the M'4 kernel coefficients
    [ax1,ax2,ax3,ax4,i1,i2,i3,i4]=interp_coeff_mp4(C(1),nx);
    % When Out of boundaries, ax might be 0
    mesh(1:nd,i1) = mesh(1:nd,i1) + ax1*part_values(1:nd,i);
    mesh(1:nd,i2) = mesh(1:nd,i2) + ax2*part_values(1:nd,i);
    mesh(1:nd,i3) = mesh(1:nd,i3) + ax3*part_values(1:nd,i);
    mesh(1:nd,i4) = mesh(1:nd,i4) + ax4*part_values(1:nd,i);
end
