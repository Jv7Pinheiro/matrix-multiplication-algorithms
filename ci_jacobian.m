function J = ci_jacobian(Mat, varargin)
%CI_JACOBIAN computes the jacobian of the problem X - Approx
%
%   J = CI_JACOBIAN(Mat) returns the jacobian of the problem X - [A, B, C],
%   where A, B, C are the composed of the cyclic invariant matrices S, U,
%   V, W in Mat.
%
%   J = CI_JACOBIAN(Mat, beta) returns the jacobian of the problem with the
%   target matrices included.

    %% Set up parameters
    params = inputParser;
    params.addParameter('beta', 0, @(x) isscalar(x) & x >= 0)
    params.parse(varargin{:});

    %% Copy from params object
    beta = params.Results.beta;


    % Extract matrices from cell
    S = Mat{1};
    U = Mat{2};
    V = Mat{3};
    W = Mat{4};
    
    % Obtain size of tensor, and cyclic invariant matrices
    n = size(S, 1);
    Rs = size(S, 2);
    Rc = size(U, 2);
    I = eye(n);

    %% Precomputations
    % Khatri-Rao Products
    SoS = khatrirao(S, S);
    VoW = khatrirao(V, W);
    WoV = khatrirao(W, V);
    WoU = khatrirao(W, U);
    UoW = khatrirao(U, W);
    UoV = khatrirao(U, V);
    VoU = khatrirao(V, U);
    
    % Permutation Matrices
    P2 = kron(eye(n),perfect_shuffle(n,n));
    P3 = perfect_shuffle(n^2,n);
    
    %% Compute Jacobian
    Js = kron(SoS, I) + P2'*kron(SoS, I) + P3'*kron(SoS, I);
    Ju = kron(VoW, I) + P2'*kron(WoV, I) + P3'*kron(VoW, I);
    Jv = kron(WoU, I) + P2'*kron(UoW, I) + P3'*kron(WoU, I);
    Jw = kron(UoV, I) + P2'*kron(VoU, I) + P3'*kron(UoV, I);
    
    % If user specifies beta, we return larger Jacobian with BI and 0 block
    % matrices below, else we return regular Jacobian
    if beta ~= 0
        BetaRows = sqrt(Beta)*eye(n*Rs + 3*n*Rc);
        J = [Js Ju Jv Jw; BetaRows];
    else
        J = [Js Ju Jv Jw];
    end
end
