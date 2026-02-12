function K = CreateTwoFixAlgorithm
% builds the algorithm with Z_3 symmetry and 2 fixed components

    % eqn (12)
    a0 = [1 0 0; 0 0 1; 0 0 1];
    
    % eqn (13)
    a1 = [0 0 0; 0 1 -1; 0 0 0];

    % eqn (14)
    U{1} = [0 0 0; 1 0 -1; 0 0 0];
    V{1} = [0 1 0; 0 0 1; 0 0 0];
    W{1} = [0 0 0; 0 1 -1; 0 1 -1];

    % eqn (15)
    U{2} = [1 0 0; 0 0 1; 0 0 0];
    V{2} = [0 1 0; 0 0 1; 0 0 1];
    W{2} = [0 0 0; 1 0 -1; 0 1 -1];

    % eqn (16)
    U{3} = [0 -1 1; 0 0 0; 0 0 0]; 
    V{3} = [0 0 0; 0 0 0; 0 0 1];
    W{3} = [0 0 0; 0 0 0; 1 0 0];

    % eqn (17)
    U{4} = [1 0 0; 1 0 0; 0 0 0];
    V{4} = [0 -1 1; 0 0 0; 0 0 0];
    W{4} = [0 0 0; 0 0 0; 0 1 0];

    % eqn (18)
    U{5} = [1 0 0; 0 0 0; 0 0 0];
    V{5} = [0 0 1; 0 0 1; 0 0 1];
    W{5} = [0 0 0; 0 0 0; 1 -1 0];

    % eqn (19)
    U{6} = [0 0 0; 0 0 1; 0 0 0];
    V{6} = [0 1 0; 0 1 0; 0 1 0];
    W{6} = [0 0 0; -1 1 0; 0 0 0];

    % eqn (20)
    U{7} = [0 0 0; 0 0 1; 0 0 1];
    V{7} = [1 0 0; 1 0 0; 1 0 0];
    W{7} = [-1 1 0; 0 0 0; 0 0 0];
    
    % symmetric part
    av = [a0(:) a1(:)];
    K = ktensor({av,av,av});
    
    % nonsymmetric part
    u = zeros(9,7); v = zeros(9,7); w = zeros(9,7);
    for i=1:7
        u(:,i) = U{i}(:);
        v(:,i) = V{i}(:);
        w(:,i) = W{i}(:);
    end
    K = K + ktensor({u,v,w});
    K = K + ktensor({w,u,v});
    K = K + ktensor({v,w,u});
    
    % col perm to obtain Z3-inv order format (S A B C)
    perm = order_for_z3inv(K.U{1},K.U{2},K.U{3});
    K = arrange(K,perm);

    % check correctness
    assert(norm(full(K) - matmul_tensor(3,3,3)) == 0);

end
