function [A,B,C] = strassen222_s1()
% generates Strassen's factor matrices for 2x2x2

    A = [
      1  0  1  0  1 -1  0 
      0  1  0  0  0  1  0 
      0  0  0  0  1  0  1 
      1  1  0  1  0  0 -1
    ];
    
    B = [
      1  1  0 -1  0  1  0 
      0  0  0  1  0  0  1 
      0  0  1  0  0  1  0 
      1  0 -1  0  1  0  1     
    ];
    
    C = [
      1  0  0  1 -1  0  1 
      0  0  1  0  1  0  0 
      0  1  0  1  0  0  0 
      1 -1  1  0  0  1  0   
    ];

    % This is the one that came from FactorStrassen()
    % U = [
    %      0     0
    %      1     0
    %      0     0
    %      1     1];
    % V = [     
    %      0     1
    %      0     0
    %      1     1
    %     -1     0];
    % W = [     
    %      1    -1
    %      0     1
    %      0     0
    %      0     0];

    % This came from the Original Strassen's Algorithm
    % U = [
    %      0     0
    %      0     0
    %      1     0
    %      1     1];
    % V = [     
    %      0     1
    %      1     1
    %      0     0
    %      -1    0];
    % W = [     
    %      1    -1
    %      0     0
    %      0     1
    %      0     0];
end