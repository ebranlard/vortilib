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

