function get_ranks(Soln, print_soln)
    nsol = length(Soln);
    nsym = size(Soln{1}{1}, 2);
    sz = size(Soln{1}{1}, 1);

    for i = 1:nsol
        rank_cnt = zeros(2, 1);
        rank_loc = zeros(nsym, 2);
        for j = 1:nsym
            mat = reshape(Soln{i}{1}(:, j), sqrt(sz), sqrt(sz));
            rnk = rank(mat);

            index = find(rank_cnt(1, :) == rnk);
            
            if ~isempty(index)
                % If the value exists, increase the count in the second row
                rank_cnt(2, index) = rank_cnt(2, index) + 1; 
            else
                % If the value doesn't exist, add a new column
                rank_cnt = [rank_cnt, [rnk; 1]];
            end

            rank_loc(j, :) = [j, rnk];
        end

        if (print_soln)
            rank_cnt(:, 2:end); % Toggle
            fprintf('\\textbf{Rs %d - Soln %d}\n', nsym, i);
            K = cyc2fac(Soln{i});
            print_char_poly(K{1}, nsym);
            % print_char_poly(Thi, nsym);
        end

        fprintf("\\begin{table}\n")
        fprintf("\\begin{tabular}{|c|c|}\n")
        fprintf("\\hline\n")
        fprintf("\\textbf{Col Number} & \\textbf{Rank}\\\\\n")
        for i = 1:nsym
            fprintf("%d & %d \\\\ \n", rank_loc(i, 1), rank_loc(i, 2))
        end
        fprintf("\\hline\n")
        fprintf("\\end{tabular}\n")
        fprintf("\\end{table}\n")
        disp('\newpage')
        fprintf("\n\n")

        % Find the column whose rank equals the square root of size

        % Find the columns whose rank is greater than 1
        
        % Perform Kronecker operation to see which of these columns commute
    end
end
