function make_matrices(K, r, s, name)

    % % use name to define build function
    % build_alg = str2func(['build_' name '_alg']);
    % % build algorithm's ktensor
    % K = build_alg();
    % 
    % % deduce matmul dimension
    dim = sqrt(size(K{1},1));
    % 
    % % compute number of symmetric components and triplets
    % r = length(K.lambda);
    % s = 0;
    % for i = 1:r
    %     if(norm(K.U{1}(:,i)-K.U{2}(:,i))+norm(K.U{2}(:,i)-K.U{3}(:,i)) == 0)
    %         s = s+1;
    %     end
    % end
    t = (r-s)/3;
    
    % open tex file to write and add header
    % fID = fopen(['../' name '_matrices.tex'],'w');
    fprintf('\\begin{align}\n\\Mthree\n&= ');
    
    % write symmetric components
    % sets number of sym comps printed on one line
    num_sym_per_line = 3;
    if strcmp(name,'twofix')
        num_sym_per_line = 1; 
    end
    for i = 1:s
        % label line
        if num_sym_per_line == 1 || mod(i,num_sym_per_line) == 1
            fprintf('\\label{%s-c%s} ', name, num2str(i));
        end      
        m = reshape(K{1}(:,i),dim,dim);
        print_matrix(m);
        fprintf('^{\\otimes 3} ');
        if i == s || mod(i,num_sym_per_line) == 0
            fprintf('\\\\\n&+ ');
        else
            fprintf('\n\t+ ');
        end
    end
    for i = 1:t
        u = K{1}(:,i+s);
        v = K{2}(:,i+s);
        w = K{3}(:,i+s);
        fprintf('\\label{%s-c%s} ', name, num2str(s+i));
        fprintf('\\mathbb{Z}_3^{std} \\cdot \n\t');
        print_matrix(reshape(K{1}(:,i+s),dim,dim));
        fprintf(' \\otimes \n\t');
        print_matrix(reshape(K{2}(:,i+s),dim,dim));
        fprintf(' \\otimes \n\t');
        print_matrix(reshape(K{3}(:,i+s),dim,dim));        
        if i == t
            fprintf(' \n');
        else
            fprintf(' \\\\ \n&+ ');
        end
    end
    
    % add footer and close tex file
    fprintf('\\end{align}\n');
    % fclose(fID);
    
    function print_matrix(m)
        dimension = size(m,1);
        fprintf('\\begin{pmatrix} ');
        for j = 1:dimension
            for k = 1:dimension
                if k < dimension
                    fprintf('%d & ', m(j,k));
                else
                    fprintf('%d ', m(j,k));
                end
            end
            if j < dimension
                fprintf('\\\\ ');
            end
        end
        fprintf('\\end{pmatrix}');
    end

end