function make_pairing(K, s, t)

    % % use name to define build function
    % build_alg = str2func(['build_' name '_alg']);
    % % build algorithm's ktensor
    % K = build_alg();
    % r = length(K.lambda);
    
    % % compute number of symmetric components and triplets
    % s = 0;
    % for i = 1:r
    %     if(norm(K.U{1}(:,i)-K.U{2}(:,i))+norm(K.U{2}(:,i)-K.U{3}(:,i)) == 0)
    %         s = s+1;
    %     end
    % end
    % t = (r-s)/3;
    
    % get rank-one component information 
    [us,vs,uapps,vapps,ui,vi] = get_rank_one_comps([K{1}, K{2}, K{3}, K{4}]);
    
    % open tex file to write and add header
    % fID = fopen(['../' name '_pairing.tex'],'w');
    fprintf('\\begin{tikzpicture}[node distance=3cm and 1cm, auto]\n\n');
    
    % print top row, corresponding to u vectors
    fprintf('\t%% top row of nodes\n');
    for i = 1:length(uapps)
        fprintf('\t\\node[main node, label={\\small %d}] (%du) ', ...
            uapps(i), i);
        if i > 1
            fprintf('[right=of %du] ',i-1);
        end
        fprintf('{$\\begin{matrix} %d & %d & %d & %d\\end{matrix}$};\n', ...
            us(1,i), us(2,i), us(3,i), us(4,i));
    end
    
    % print bottom row, corresponding to v vectors
    fprintf('\n\t%% bottom row of nodes\n');
    for j = 1:length(vapps)
        fprintf('\t\\node[main node, label=below:{\\small %d}] (%dv) [below=of %du] ', ...
            vapps(j), j, j);
        fprintf('{$\\begin{matrix} %d & %d & %d & %d \\end{matrix}$};\n', ...
            vs(1,j), vs(2,j), vs(3,j), vs(4,j));
    end
    
    % define edge colors and ways to bend for multi-edges
    colors = {'black','blue','red','green','olive','purple','orange','cyan','brown', ...
        'blue!50','red!50','green!50','olive!50','purple!50','orange!50','cyan!50','brown!50'};
    bends = {'','bend left=8','bend right=8','bend left=16','bend right=16'};
    edgecount = ones(size(us,2),size(vs,2));
    
    % print edges
    fprintf('\n\t%% edges between rows\n\t\\path[line]\n');
    % loop over symmetric components
    for i = 1:s
        % check if matrix is rank one
        if ui(i) && vi(i)
            fprintf('\t(%du) edge [dashed,%s,%s] node {} (%dv)\n',ui(i),colors{1},bends{edgecount(ui(i),vi(i))},vi(i));
            edgecount(ui(i),vi(i)) = edgecount(ui(i),vi(i)) + 1;
        end
    end
    % loop over triplets
    for i = 1:t
        % loop within triplet
        for j = 0:2
            % check if matrix is rank one
            if ui(i+s+j*t) && vi(i+s+j*t)
                fprintf('\t(%du) edge [%s,%s] node {} (%dv)\n',ui(i+s+j*t),colors{1+i},bends{edgecount(ui(i+s+j*t),vi(i+s+j*t))},vi(i+s+j*t));
                edgecount(ui(i+s+j*t),vi(i+s+j*t)) = edgecount(ui(i+s+j*t),vi(i+s+j*t)) + 1;
            end
        end
    end
    fprintf('\t;\n');
    
    % add footer and close tex file
    fprintf('\n\\end{tikzpicture}\n');
    % fclose(fID);
end