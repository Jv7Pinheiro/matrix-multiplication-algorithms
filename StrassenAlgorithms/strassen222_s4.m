function [U,V,W] = strassen222_s4() 
% return U,V,W representation of 2x2 rank 7 decomposition with 4 sym
% from Chiantini, Landsberg et al 2016 paper
% U and V are col-major, W is row-major

    S = zeros(4,4);
    A = zeros(4,1);
    B = zeros(4,1);
    C = zeros(4,1);
    
    S(2,1)=-1; S(3,1)=1; S(4,1)=-1;
    S(1,2)=1;
    S(2,3)=1; S(4,3)=1;
    S(3,4)=-1; S(4,4)=1;

    A(2,1)=1;
    B(3,1)=1;
    C=[1 -1 1 -1]';

    U = [S A B C];
    V = [S C A B];
    W = [S B C A];
    
end