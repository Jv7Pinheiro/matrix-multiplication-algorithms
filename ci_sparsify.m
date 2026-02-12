function [Q, innnz, outnnz] = ci_sparsify(P, threshold, varargin)
%CI_SPARSIFY introduces zeros while maintaining the same error
%
%   Q = CI_SPARSIFY(P, threshold) uses the threshold to define how many
%   entries are considered nonzeros, which are entries less than the
%   threshold. It then selects the best column in each of the cyclic
%   matrices to introduce zeros to each cylic matrix respectively.  
%
%   Q = CI_SPARSIFY(P, threshold, 'param', value, ...) pecifies additional
%   parameters for the function. Specifically...
%   
%   'print' - prints information regardning errors and number of nonzeros
%   'type' - chooses the method that will introduce the nonzeros, which can
%   be one of two of the following strings:
%       'schur'  schur decomposition
%       'eig'    eigenvalue decomposition

    %% Set up parameters
    params = inputParser;
    params.addParameter('type', 'schur', @(x) isstring(x))
    params.addParameter('print', 0, @isscalar);
    params.parse(varargin{:});

    %% Copy from params object
    type = params.Results.type;
    print = params.Results.print;

    type = lower(type);
    if ~(isequal(type, 'eig') || isequal(type, 'schur'))
        error('function type must be either eig or schur');
    end

    %Q = P;
    dim = sqrt(size(P{1},1));
    Rs = size(P{1},2);
    Rc = size(P{2},2);
    
    %tmp
    T = matmul_tensor(dim,dim,dim);
    
    % check input error
    U = cyc2fac(P);
    inerr = norm(T-full(ktensor(U)))^2;
    
    % count input nnz
    if Rs == 0
        innnz = sum(cellfun(@(x) nnz(abs(x) > threshold), {P{2}, P{3}, P{4}}));
    else
        innnz = sum(cellfun(@(x) nnz(abs(x) > threshold), P));
    end
    
    if (print == 1)
        fprintf('\nInitial error is %7.1e, with %d zeros\n', inerr, innnz);
    end

    % try symmetric components
    if Rs ~= 0
        Qs = cell(Rs,1);
        nnzs = zeros(Rs,1);
        if isequal(type, 'schur')
            for i = 1:Rs
                [Qs{i},~] = schur(reshape(P{1}(:,i),dim,dim));
                Q = cellfun(@(x) kron(Qs{i}',Qs{i}') * x, P, 'UniformOutput', false);
                nnzs(i) = sum(cellfun(@(x) nnz(abs(x) > threshold), Q)); 
            end
        elseif (isequal(type, 'eig'))
            for i = 1:Rs
                PP = (reshape(P{1}(:, i),dim,dim)' + reshape(P{1}(:, i),dim,dim))/2;
                [Qs{i},~] = eig(reshape(PP,dim,dim));
                Q = cellfun(@(x) kron(Qs{i}',Qs{i}') * x, P, 'UniformOutput', false);
                nnzs(i) = sum(cellfun(@(x) nnz(abs(x) > threshold), Q)); 
            end
        end
    end

    % try cyclic components
    if Rc ~= 0
        Qc = cell(Rc,3);
        nnzc = zeros(Rc,3);
        if isequal(type,'schur')
            for i = 1:Rc
                for j = 1:3
                    [Qc{i,j},~] = schur(reshape(P{j+1}(:,i),dim,dim));
                    Q = cellfun(@(x) kron(Qc{i,j}',Qc{i,j}') * x, P, 'UniformOutput', false);
                    nnzc(i,j) = sum(cellfun(@(x) nnz(abs(x) > threshold), Q));
                end
            end
        elseif isequal(type, 'eig')
            for i = 1:Rc
                for j = 1:3
                    PP = (reshape(P{j+1}(:, i),dim,dim)' + reshape(P{j+1}(:, i),dim,dim))/2;
                    [Qc{i,j},~] = eig(reshape(PP,dim,dim));
                    Q = cellfun(@(x) kron(Qc{i,j}',Qc{i,j}') * x, P, 'UniformOutput', false);
                    nnzc(i,j) = sum(cellfun(@(x) nnz(abs(x) > threshold), Q));
                end
            end
        end
    end
    
    % choose sparsest solutiom
    if Rs == 0
        min_nnzc = min(nnzc(:));
        [i_nnzc,j_nnzc] = find(nnzc == min_nnzc,1);
        Q = cellfun(@(x) kron(Qc{i_nnzc,j_nnzc}', Qc{i_nnzc,j_nnzc}') * x, P, 'UniformOutput', false);
    else 
        [min_nnzs,i_nnzs] = min(nnzs);
        min_nnzc = min(nnzc(:));
        if min_nnzs < min_nnzc
            Q = cellfun(@(x) kron(Qs{i_nnzs}',Qs{i_nnzs}') * x, P, 'UniformOutput', false);
        else
            [i_nnzc,j_nnzc] = find(nnzc == min_nnzc,1);
            Q = cellfun(@(x) kron(Qc{i_nnzc,j_nnzc}', Qc{i_nnzc,j_nnzc}') * x, P, 'UniformOutput', false);
        end
    end

    % check output error
    Unew = cyc2fac(Q);
    outerr = norm(T-(ktensor(Unew)))^2;

    % count output nnz
    if Rs ~= 0
        outnnz = sum(cellfun(@(x) nnz(abs(x) > threshold), Q));
    else
        outnnz = sum(cellfun(@(x) nnz(abs(x) > threshold), {Q{2}, Q{3}, Q{4}}));
    end

    if (print == 1)
        fprintf('Output error is %7.1e, with %d zeros\n\n', outerr, outnnz);
    end

end