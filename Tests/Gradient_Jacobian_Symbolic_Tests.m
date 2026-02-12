%% Gradient
addpath("C:\Users\jv777\OneDrive\Documentos\Documents\Wake Forest University - Graduate School Research\tensor_toolbox-v3.5")
% Defining the model using symbolic vars
n = 4;
Rs = 1;
Rc = 2;
Rank = Rs + 3*Rc;

Beta = sym('Beta', [1, 1], 'int');

S = sym('S', [n, Rs], 'int');
U = sym('U',[n, Rc], 'int');
V = sym('V',[n, Rc], 'int');
W = sym('W',[n, Rc], 'int');

St = sym('St', [n, Rs], 'int');
Ut = sym('Ut',[n, Rc], 'int');
Vt = sym('Vt',[n, Rc], 'int');
Wt = sym('Wt',[n, Rc], 'int');

X = sym('X',[n, n, n], 'int');

A = [S U V W];
B = [S W U V];
C = [S V W U];

% sum of outer products column vectors
X_hat = zeros(n, n, n);

for i=1:Rank
    X_hat = X_hat + reshape(kron(kron(C(:,i), B(:,i)), A(:,i)), n, n, n);
end

phi = [X(:)-X_hat(:); Beta*(S(:)-St(:)); Beta*(U(:)-Ut(:)); Beta*(V(:)-Vt(:)); Beta*(W(:)-Wt(:))]; %phi
r = sum((phi .*phi)/2, "all");

%Compare Gradients
for i=1:1
    % randomize values to test gradient differences
    S_vals = randi(10, n, Rs);
    U_vals = randi(10, n, Rc);
    V_vals = randi(10, n, Rc);
    W_vals = randi(10, n, Rc);

    St_vals = randi(10, n, Rs);
    Ut_vals = randi(10, n, Rc);
    Vt_vals = randi(10, n, Rc);
    Wt_vals = randi(10, n, Rc);

    X_vals = randi(10, n, n, n);

    % put into factor matrices
    FM = {S_vals, U_vals, V_vals, W_vals};
    FtM = {St_vals, Ut_vals, Vt_vals, Wt_vals};

    %Just pass value matrices, and compute four gradients

    %compute gradient 1 from function
    [gradient1, f] = ci_gradient_functionvalue(FM, FtM, tensor(X_vals), 9);

    % compute gradient 2 based on sym function and subbing variables
    gradient_sym = gradient(r, [S(:); U(:); V(:); W(:)]);
    gradient2 = subs(gradient_sym, S, S_vals);
    gradient2 = subs(gradient2, St, St_vals);
    gradient2 = subs(gradient2, U, U_vals);
    gradient2 = subs(gradient2, Ut, Ut_vals);
    gradient2 = subs(gradient2, V, V_vals);
    gradient2 = subs(gradient2, Vt, Vt_vals);
    gradient2 = subs(gradient2, W, W_vals);
    gradient2 = subs(gradient2, Wt, Wt_vals);
    gradient2 = subs(gradient2, X, X_vals);
    gradient2 = subs(gradient2, Beta, 3);
    gradient2 = cast(gradient2, "double");

    S_symval = reshape(gradient2(1:n*Rs),n,Rs);
    U_symval = reshape(gradient2(n*Rs+1:n*Rs+n*Rc),n,Rc);
    V_symval = reshape(gradient2(n*Rs+1 + n*Rc:n*Rs+2*n*Rc),n,Rc);
    W_symval = reshape(gradient2(n*Rs+1 +2*n*Rc:n*Rs+3*n*Rc),n,Rc);


    % message if gradients are ever differnt
    assert(isequal(gradient1, gradient2), "Gradients are different\n")
end
fprintf("Gradients are the same\n")


%% Jacobian
%Compare Jacobians
for i=1:100
    % randomize values to hessians
    S_vals = randi(10, n, Rs);
    U_vals = randi(10, n, Rc);
    V_vals = randi(10, n, Rc);
    W_vals = randi(10, n, Rc);
    X_vals = randi(10, n, n, n);

    % put into factor matrices
    FM = {S_vals, U_vals, V_vals, W_vals};

    %compute jacobian1 from our minres function
    jacobian1 = ci_jacobian(FM, 9);


    % compute jacobian 2 based on sym functions and subbing variables
    jacobian_sym = jacobian(phi(:), [S(:); U(:); V(:); W(:)]);

    jacobian2 = subs(jacobian_sym, S, S_vals);
    jacobian2 = subs(jacobian2, U, U_vals);
    jacobian2 = subs(jacobian2, V, V_vals);
    jacobian2 = subs(jacobian2, W, W_vals);
    jacobian2 = subs(jacobian2, Beta, -3);
    jacobian2 = subs(jacobian2, X, X_vals);
    jacobian2 = cast(jacobian2, "double");
    
    diff = jacobian1 + jacobian2;
    % message if jacobians are ever differnt
    assert(isequal(jacobian1, -jacobian2), "Jacobians are different")
end
fprintf("Jacobians are the same\n")