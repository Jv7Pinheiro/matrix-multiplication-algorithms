function [U,V,W] = strassen444_s1()
% return U,V,W representation of 4x4 rank 49 decomposition with 16 sym
% from Kronecker'ing Strassen soln for 2x2
% U and V are col-major, W is row-major

    % cyclic invariant structure of Strassen
    S = [1 0 0 1]';
    A = [
         0     0
         1     0
         0     0
         1     1];
    B = [
         0     1
         0     0
         1     1
        -1     0];
    C = [
         1    -1
         0     1
         0     0
         0     0];

    % create permutation matrix that swaps Morton for col-ordering
    P = eye(16);
    P(3,3)=0;P(5,5)=0;P(3,5)=1;P(5,3)=1;
    P(4,4)=0;P(6,6)=0;P(4,6)=1;P(6,4)=1;
    P(11,11)=0;P(13,13)=0;P(11,13)=1;P(13,11)=1;
    P(12,12)=0;P(14,14)=0;P(12,14)=1;P(14,12)=1;

    % build U V W in cyclic-invariant structure
    S4 = kron(S,S);
    A4 = [kron(S,A) kron(A,S) kron(A,A) kron(A,B) kron(A,C)];
    C4 = [kron(S,C) kron(C,S) kron(C,C) kron(C,A) kron(C,B)];
    B4 = [kron(S,B) kron(B,S) kron(B,B) kron(B,C) kron(B,A)];
    U = P*[S4 A4 B4 C4];
    V = P*[S4 C4 A4 B4];
    W = P*[S4 B4 C4 A4];
end