function A = vec2cyc(x, Shape)
%VEC2CYC Reshapes the vector x into cyclic matrices S U V W
%   
%   A = vec2cyc(x, SHAPE) reshapes the vector x into the factor matrices
%   whose sizes are taken from the cell array SHAPE

    n_sqr = size(Shape{1}, 1);
    Rs = size(Shape{1}, 2);
    Rc = size(Shape{2}, 2);
    
    A = cell(4,1);
    for i = 1:4
        if (i == 1)
            A{i} = reshape(x(1:n_sqr*Rs), [], Rs);
        else
            idx1 = n_sqr*Rs + (i-2)*n_sqr*Rc + 1;
            idx2 = n_sqr*Rs + (i-1)*n_sqr*Rc;
            A{i} = reshape(x(idx1:idx2), [], Rc);
        end
    end

end

