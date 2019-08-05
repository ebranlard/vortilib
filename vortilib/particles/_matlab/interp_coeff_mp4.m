function [ax1,ax2,ax3,ax4,i1,i2,i3,i4]=interp_coeff_mp4(i2,dx2,nx)
    if nargin==0
        % TEST
        [ax1,ax2,ax3,ax4,i1,i2,i3,i4]=interp_coeff_mp4(3,0.1,5)
        return
    end
    % !---------------------------------------------------------------------------
    % ! Find index of other cells
    % !---------------------------------------------------------------------------
    i1 = i2 - 1;
    i3 = i2 + 1;
    i4 = i2 + 2;
    % !---------------------------------------------------------------------------
    % ! Find normalised distance to other cells
    % !---------------------------------------------------------------------------
    dx1 = dx2 + 1.0;
    dx3 = 1.0 - dx2 ;
    dx4 = 2.0 - dx2 ;
    % !---------------------------------------------------------------------------
    % ! The M-prime-4 kernel
    % !---------------------------------------------------------------------------
    ax1 = 0.5 * (2.0 - dx1)^2 * (1.0 - dx1);
    ax2 = 1.0 - 2.5*dx2^2 + 1.5*dx2^3;
    ax3 = 1.0 - 2.5*dx3^2 + 1.5*dx3^3;
    ax4 = 0.5 * (2.0 - dx4)^2 * (1.0 - dx4);

    if(i2==nx-1) 
        % My hack for one sided right
%         ax3 = 2*(ax3)+ax4;
%         ax2 = (ax2+ax4);
%         i4=1; ax4=0;% arbitrary, i1 will not be used
        i1=1; i2=1; i3=1; i4=1;
        ax1=0; ax2=0; ax3=0; ax4=0;
    elseif (i2==1) 
%         % My hack for one sided left
%         ax2 = 2*(ax2)+ax1;
%         ax3 = (ax3+ax1);
%         i1=1; ax1=0; % arbitrary, i1 will not be used
        i1=1; i2=1; i3=1; i4=1;
        ax1=0; ax2=0; ax3=0; ax4=0;
    elseif (i2<1) 
        % Should not happen
        i1=1; i2=1; i3=1; i4=1;
        ax1=0; ax2=0; ax3=0; ax4=0;
    elseif (i2>nx-1) 
        % Should not happen
        i1=1; i2=1; i3=1; i4=1;
        ax1=0; ax2=0; ax3=0; ax4=0;
    elseif (i1<=0 || i2<=0) 
        % Might happen if on grid
        i1=1; i2=1; i3=1; i4=1;
        ax1=0; ax2=0; ax3=0; ax4=0;
    elseif (i4>nx || i3>nx) 
        % Might happen if on grid
        i1=1; i2=1; i3=1; i4=1;
        ax1=0; ax2=0; ax3=0; ax4=0;
    end

end
