function [U, V, W] = strassen444_s4_v2()
% return rank-49 soln based on the kronecker product of
% 222 alg with 1 sym comp and 222 alg with 4 sym comps
% U and V are col-major, W is row-major

    [u1,v1,w1] = strassen222_s4;
    [u2,v2,w2] = strassen222_s1;

    % create permutation matrix that swaps Morton for col-ordering
    P = eye(16);
    P(3,3)=0;P(5,5)=0;P(3,5)=1;P(5,3)=1;
    P(4,4)=0;P(6,6)=0;P(4,6)=1;P(6,4)=1;
    P(11,11)=0;P(13,13)=0;P(11,13)=1;P(13,11)=1;
    P(12,12)=0;P(14,14)=0;P(12,14)=1;P(14,12)=1;

    % create permutation matrix that reorders symmetric components
    % to be first followed by blocks of triples
    Pc = eye(49);
    order=[1,8,15,22, ...       % s
            2,4,9,11,16,18,23,25,29,30,31,32,33,34,35, ...  % a
            7,5,14,12,21,19,28,26,36,42,37,40,41,39,38, ...  % b
            3,6,10,13,17,20,24,27,43,45,49,48,46,47,44  ];     % c
    Pc=Pc(:,order);
    
    U = P*kron(u1,u2)*Pc;
    V = P*kron(v1,v2)*Pc;
    W = P*kron(w1,w2)*Pc;
    
end