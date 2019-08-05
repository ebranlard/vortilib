function [ic dc]=fCoordRegularGrid(x0,xBox,dBox)
    if nargin==0
        [ic,dc]=fCoordRegularGrid(1.1,0,2)
        return
    end
    C = (x0-xBox)/dBox; % position in grid coordinates
    ic  = floor(C)+1 ;  % front lower left grid point for the cell containing the point
    dc = C+1-ic     ;  % Complementary Distance to the front lower left node
