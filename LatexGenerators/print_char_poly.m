function P = print_char_poly(U,s,varargin)
% assumes Z3-invariant solutions, so U's are enough info
% s is number of symmetric components, 
% optional 3rd argument is name for printing
% returns characteristic polynomials P and
% prints latex table of char poly representations

    if nargin < 3
        name = 'name';
    else
        name = varargin{1};
    end
    
    % make sure U is a cell array
    if ~iscell(U)
        tmp{1} = U;
        U=tmp;
    end

    % deduce matmul dim (assume square) and rank
    dim = round(sqrt(size(U{1},1)));
    rnk = size(U{1},2);
    triples = (rnk-s)/3;
    
    % initialize eigenvalues matrix, char polynomial arrays, and symbolic variable
    P = sym('p',[s,1]);
    Q = sym('q',[triples,3]);
    syms t;
    
    % table header
    % disp('\begin{table}');
    disp('\begin{tabular}{|c|c|c|}');
    disp('\hline');
    disp('\textbf{Type} & \textbf{Char.~Poly.} & \textbf{Count} \\');
    
    
    % loop over algorithms
    for a = 1:1 % length(U)
        
        % algorithm info
        fprintf('\\hline \\multicolumn{3}{|c|}{\\texttt{%s}} \\\\ \\hline\n', name);
       
        % loop over symmetric components (assumed to be first s columns)
        for i = 1:s
            M = reshape(U{a}(:,i),dim,dim);
            % get characteristic polynomial
            P(i) = simplify(det(t*eye(dim)-M));
        end  

        [C,~,IC] = unique(P,'rows');
        nu = size(C,1);
        fprintf('\\multirow{%s}{*}{symmetric}\n', num2str(nu));
        for j = 1:nu
            pp = strrep(char(C(j)),'*','');
            fprintf('& $%s$ & %s \\\\\n',pp,num2str(nnz(IC==j)));
        end
        
        disp('\hline');
                
        % loop over non-symmetric components
        for i = 1:triples
            for k = 1:3
                M = reshape(U{a}(:,s+i+(k-1)*triples),dim,dim);
                % get characteristic polynomial
                Q(i,k) = simplify(det(t*eye(dim)-M));
            end
        end 
        
        % impose ordering on triples by sorting rows of Q
        Q = sort(Q,2);

        [C,~,IC] = unique(Q,'rows');
        nu = size(C,1);
        disp(['\multirow{',num2str(nu),'}{*}{triples}']);
        for j = 1:nu
            disp([' & $\left\{',strrep(char(C(j,1)),'*',''),',',strrep(char(C(j,2)),'*',''),',',strrep(char(C(j,3)),'*',''),'\right\}$ & ',num2str(nnz(IC==j)),' \\']);
        end
    
    end
    
    % table footer
    disp('\hline');
    disp('\end{tabular}');
    % disp('\end{table}');
    
end