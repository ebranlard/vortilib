% test cases
function vectoranalysis()
% NOTE:  DO NOT CHANGE OR CHANGE vectoranalysis.py accordingly
clc;


% --- Gradient  1d
% x = [0,1,3,4];
% x2 = x.^2;
% fx = gradient(x2)
% fx = gradient(x2,2)
% fx = gradient(x2,x)
% fx = gradient(x2')
% fx = gradient(x2',x)
% fx = gradient(x2',x')

% --- Gradient  2d
% x = -1:0.5:1.5
% y = (0:0.2:1.2)'
% [X,Y]=meshgrid(x,y)
% Z = X .* exp(-X.^2 - Y.^2)
% 
% [px,py] = gradient(Z);
% pythondisp(px(1:2,1:2),'px');
% pythondisp(py(1:2,1:2),'py');
% 
% [px,py] = gradient(Z,0.5);
% pythondisp(px(1:2,1:2),'px');
% pythondisp(py(1:2,1:2),'py');
% 
% [px,py] = gradient(Z,0.5,0.2);
% pythondisp(px(1:2,1:2),'px');
% pythondisp(py(1:2,1:2),'py');
% 
% [px,py] = gradient(Z,x,y);
% pythondisp(px(1:2,1:2),'px_ref');
% pythondisp(py(1:2,1:2),'py_ref');


% --- Divergence 2d
% x = -1:0.5:1.5;
% y = (0:0.2:1.2)';
% [X,Y]=meshgrid(x,y);
% U =    X .* exp(-X.^2 - Y.^2)
% V = Y.^2 .* exp(-X.^2 - Y.^2)
% 
% d=divergence(U,V);
% pythondisp(d(1:2,1:2),'div_ref');
% 
% d=divergence(X,Y,U,V)
% pythondisp(d(1:2,1:2),'div_ref');

% --- Curl 2d
x = -1:0.5:1.5;
y = (0:0.2:1.2)';
[X,Y]=meshgrid(x,y);
U =    X .* exp(-X.^2 - Y.^2)
V = Y.^2 .* exp(-X.^2 - Y.^2)

[c,av]=curl(U,V);
pythondisp(c(1:2,1:2),'curl_ref');

[c,av]=curl(X,Y,U,V);
pythondisp(c(1:2,1:2),'curl_ref');
c




end


function pythondisp(m,v)
    fprintf('%s = np.array([\n',v)
    for i=1:size(m,1)
        fprintf('[')
        for j =1:size(m,2)
            fprintf('%12.7f',m(i,j)) 
            if j<size(m,2)
                fprintf(', ')
            end
        end
        fprintf(']')
        if i<size(m,1)
            fprintf(',\n')
        end
    end
    fprintf('])\n')
end

