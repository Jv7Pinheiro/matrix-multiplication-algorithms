function P = PerfectShuffle(m,n)
%PERFECT_SHUFFLE returns perfect shuffle permutation matrix that reverses 
%col-major and row-major vec orderings so that P*vec(A) = vec(A') for mxn 
%matrix A
%
%   P = PERFECT_SHUFFLE(m, n)

    P = zeros(m*n,m*n);
    for i = 1:m
        for j = 1:n
            P((i-1)*n+j,i+(j-1)*m) = 1;
        end
    end

end