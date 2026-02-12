function [us,vs,uapps,vapps,ui,vi] = GetRankOneComponents(U)
% assumes Z3-invariant solution, so U is factor matrix, s is # sym comps
% returns ranks of components, and 
%  if rank one, gives outer-product vectors in us and vs
%    us/vs stores vectors as columns
%    uapps/vapps stores count of appearances of each u/v
%    ui/vi stores index of u and v for each rank-one matrix (0 otherwise)

    % deduce matmul dim (assume square) and rank
    dim = round(sqrt(size(U,1)));
    rnk = size(U,2);
    
    % initialize ranks matrix and rank-one terms' vectors
    ranks = zeros(1,rnk);
    uu = zeros(dim,rnk); 
    vv = zeros(dim,rnk);
    ui = zeros(1,rnk);
    vi = zeros(1,rnk);
    
    % compute rank-one factorizations
    for i = 1:rnk
        % reshape each component
        Um = reshape(U(:,i),dim,dim);
        % compute rank
        ranks(i) = rank(Um);
        if ranks(i) == 1
            % find unit-norm vectors
            [u,~,v] = svd(Um);
            % scale to be {0,1,-1} vectors
            uu(:,i) = round( u(:,1) / max(abs( u(:,1) )) );
            vv(:,i) = round( v(:,1)  / max(abs( v(:,1) )) );
        end
    end  
    
    % get unique us (up to sign)
    us = zeros(dim,0);
    uapps = zeros(1,0);
    for i=1:rnk
        nu = size(us,2);
        % if nonzero vector
        if max(abs(uu(:,i)-zeros(dim,1)))
            newu = 1;
            for j=1:nu
                % disregard if matches existing vector (or its negative)
                if ~max(abs(uu(:,i)-us(:,j))) || ~max(abs(uu(:,i)+us(:,j)))
                    newu = 0;
                    uapps(j) = uapps(j)+1;
                end
            end
            if nu == 0 || newu == 1
                us(:,nu+1) = uu(:,i);
                uapps(nu+1) = 1;
            end
        end
    end
    
    % get unique vs (up to sign)
    vs = zeros(dim,0);
    vapps = zeros(1,0);
    for i=1:rnk
        nv = size(vs,2);
        % if nonzero vector
        if max(abs(vv(:,i)-zeros(dim,1)))
            newv = 1;
            for j=1:nv
                % disregard if matches existing vector (or its negative)
                if ~max(abs(vv(:,i)-vs(:,j))) || ~max(abs(vv(:,i)+vs(:,j)))
                    newv = 0;
                    vapps(j) = vapps(j)+1;
                end
            end
            if nv == 0 || newv == 1
                vs(:,nv+1) = vv(:,i);
                vapps(nv+1) = 1;
            end
        end
    end
    
    % flip signs to be more natural
    for i=1:length(uapps)
        % more pos than neg
        if nnz(us(:,i)<0) > nnz(us(:,i)>0)
            us(:,i) = -us(:,i);
        end
        % start with pos
        if nnz(us(:,i)) == 2 && nnz(us(:,i)<0) == 1
            if us(1,i) == -1 || (us(1,i) == 0 && us(2,i) == -1)
                us(:,i) = -us(:,i);
            end
        end
    end
    for i=1:length(vapps)
        % more pos than neg
        if nnz(vs(:,i)<0) > nnz(vs(:,i)>0)
            vs(:,i) = -vs(:,i);
        end
        % start with pos
        if nnz(vs(:,i)) == 2 &&  nnz(vs(:,i)<0) == 1
            if vs(1,i) == -1 || (vs(1,i) == 0 && vs(2,i) == -1)
                vs(:,i) = -vs(:,i);
            end
        end
    end
    
    % sort us and vs in decreasing order of appearances, then
    %  increasing number of nonzeros, then by values
    % get nnzs per column
    unnzs = ones(1,dim)*abs(us);
    vnnzs = ones(1,dim)*abs(vs);
    
    % sort by 4th index (apps), then 5th (nnzs), then values (1st,2nd,3rd)
    [~,uorder] = sortrows([us; uapps; unnzs]',[-(dim+1),dim+2,-(1:dim)]);
    [~,vorder] = sortrows([vs; vapps; vnnzs]',[-(dim+1),dim+2,-(1:dim)]);
    
    % reorder everything
    us = us(:,uorder);
    uapps = uapps(uorder);
    vs = vs(:,vorder);
    vapps = vapps(vorder);

    % loop over U components to compute ui/vi
    for i=1:rnk
        if ranks(i) == 1
            for j=1:size(us,2)
                if ~max(abs(uu(:,i)-us(:,j))) || ~max(abs(uu(:,i)+us(:,j)))
                    ui(i) = j;
                end
            end
            for j=1:size(vs,2)
                if ~max(abs(vv(:,i)-vs(:,j))) || ~max(abs(vv(:,i)+vs(:,j)))
                    vi(i) = j;
                end
            end
        end
    end
    
end

    