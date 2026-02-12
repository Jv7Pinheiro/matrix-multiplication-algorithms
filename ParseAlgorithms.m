%% Set Up
clear
clc

Rank = 49;
MatMul = 4;
Rs_list = [16]; % (Rs 1, Rc 16), (Rs 4, Rc 15), (Rs 13, Rc 12), (Rs 16, Rc 11)
% The following Rs 16 algorithms have the special property that all cyclic components map to something else
% 10 11 14 19 21 23 25 28 31


%% Loop chosen algorithm types
for Rs = Rs_list
    Rc = (Rank - Rs)/3;
    
    % Some of the .mat files are the original Strassen algorithms
    % The following boolean is a toggle in case we want to investigate them
    is_Strassen = false;
    if is_Strassen == true
        file = "Rs" + Rs + "_Rc" + Rc + "_Strassen";
    else
        file = "Rs" + Rs + "_Rc" + Rc + "_Solutions";
    end
    
    % Load algorithm
    Algs = load(file + ".mat");
    Algs = Algs.(file);
    
    % Go through each different algorithm we have for a given Rs value
    fprintf("### Rs = %d ###\n", Rs);
    
    % for i = 1:numel(Algs) % Loop through all algorithms
    % for i = [10 11 14 19 21 23 25 28 31] % Loop through selected algorithms
    for i = 10:10 % Investigate a single algorithms'
        %% Obtain information regarding symmetric components
        fprintf("\t## Algorithm %d ##\n", i);
        
        % Rs_Mats is a three way tensor where the ith frontal slice is the
        % ith component of the symmetric component reshaped into a matrix
        Rs_Mats = zeros(MatMul, MatMul, Rs); 

        % Rs_char_polys is the array containing the characteristic
        % polynomial of each symmetric component (reshaped as a matrix)
        Rs_char_polys = zeros(Rs, MatMul+1);

        % Loop through the symmetric components obtaining each char poly
        for j = 1:Rs
            Rs_Mats(:, :, j) = reshape(Algs{i}{1}(:, j), MatMul, MatMul);
            Rs_char_polys(j, :) = charpoly(Rs_Mats(:, :, j));
        end
        
        Rs_char_polys; % Print symmetric characteristic polynomials
    
        % Given the char polys of all symmetric components,
        % find all of the distinct character polynomials
        [un_cp, ~, I] = unique(Rs_char_polys, "rows");
    
        % Get the count of how many times each char poly shows up
        a_counts = accumarray(I,1);
        
        % Get the index of all of the unique cp that appear only once
        idx = find(a_counts == 1);
        Rs_unique_idx = zeros(1, length(idx));

        % In case there is more than one singleton
        for j = 1:length(idx)
            Rs_unique_idx(j) = find(all(Rs_char_polys == un_cp(idx(j), :), 2));
        end
    
        Rs_unique_idx;

        %% This section is concerned with fixing the transformation of 
        % symmetric components that must map to themselves
        
        % Rs_Maps is an Rs by Rs matrix where the (i, j) entry is a 
        %   0 if i cannot map to j
        %   1 if i can map to j
        %   2 is i can ONLY map to i
        Rs_Maps = zeros(Rs);

        % the null space of stacked_fixed_sc contains the solution to the 
        % problem I(x)A - A'(x)I = 0 where A is the Rs_Mats of the 
        % symmetric component with unique characteristic polynomial
        stacked_fixed_sc = [];

        % We loop through the unique symmetric components stacking each 
        % I(x)A - A'(x)I matrix
        for j = Rs_unique_idx
            A = Rs_Mats(:, :, j);
            Rs_Maps(j, j) = 2;
            C = kron(eye(MatMul), A) - kron(A', eye(MatMul));
            stacked_fixed_sc = [stacked_fixed_sc; C]; % Append the computed matrix C to stacked_A
        end
        
        % The rank of the null space is the number of 0 eigenvalues
        [U, S, V] = svd(stacked_fixed_sc);
        diag(S);
    
        % there will always be one as the identity commutes with itself,
        % so we ignore this case and move on to next algorithm
        if (nnz(diag(S) > 1e-13) == (MatMul^2 - 1))
            fprintf("\t\tAlgorithm does not have a diagonal action because the fixed points don't have a shared commutator\n");
            continue
        end
        
        %% This section investigates which other symmetric components map 
        % to what given that the unique ones must map to themselves

        % Loop through symmetric components veryfying the possible maps of
        % each one
        for p = 1:Rs
            if ismember(p, Rs_unique_idx) % Ignore the unique components
                continue
            end
            
            % Only check mapping to other components that have the same cp
            same_cp_idx = find(all(Rs_char_polys == Rs_char_polys(p, :), 2));
            for q = same_cp_idx'
                Bp = Rs_Mats(:, :, p);
                Bq = Rs_Mats(:, :, q);
                D = [stacked_fixed_sc; kron(eye(MatMul), Bp) - kron(Bq', eye(MatMul))];
                
                [U, S, V] = svd(D);

                % if (p == q)
                % [p, q]
                % diag(S)
                % V(:, end-3:end)
                % end
                
                % If the current check is self maps, then we must ignore
                % the identity that is always there
                if (p == q && nnz(diag(S) > 1e-13) < (MatMul^2 - 1))
                    Rs_Maps(p, q) = 1;
                end

                % if current check are two distinct components, then any
                % map is enough
                if (p ~= q && nnz(diag(S) > 1e-13) < (MatMul^2))
                    Rs_Maps(p, q) = 1;
                end
            end
        end
        
        %% This section prints information regarding the maps of the symmetric components
        
        % Print Rs_Maps
        Rs_Maps;

        % Check for symmetry
        fprintf("\t\tSymmetric Components (RS)\n")
        if (issymmetric(Rs_Maps))
            fprintf("\t\t\tRs Maps are symmetric\n")
        else
            fprintf("\t\t\tRs Maps are NOT symmetric\n")
        end
        
        % Check whether all components map to themselves
        all_self_map = false;
        if (all(diag(Rs_Maps) >= 1))
            fprintf("\t\t\tAll symmetric components map to themselves\n")
            all_self_map = true;
        else
            fprintf("\t\t\tNot all symmetric components map to themselves\n")
        end
        
        % Check whether a component maps to absolutely nothing
        % This is a stronger statement than the next check below
        has_zero_row = any(all(Rs_Maps == 0, 2));
        if (has_zero_row)
            fprintf("\t\t\tA symmetric component maps to nothing\n")
        else
            fprintf("\t\t\tAll symmetric components maps to something\n")
        end
        
        % Check whether a component doesn't map to any other component
        % This is different from above, the above check if there are any
        % zero rows (i.e. a component maps to nothing), this one removes
        % the diagonal and checks for any zero row (i.e. a component that
        % maps to no other component)
        mask = ~eye(Rs);
        off_diag = reshape(Rs_Maps(mask), Rs-1, Rs)';
        off_diag(Rs_unique_idx, :) = [];
        has_zero_row_off_diag = any(all(off_diag == 0, 2));
        if (has_zero_row_off_diag)
            fprintf("\t\t\tA symmetric component maps to no other component\n")
        else
            fprintf("\t\t\tAll symmetric components maps to some other component\n")
        end
    
        fprintf("\t\t    NNZ of Rs Maps: %d\n", nnz(Rs_Maps));

        % An idea for later on is that if every component maps to
        % themselves, then we fix the ones that only map to themselves,
        % this is different from where Rs_Maps = 2 as those only map to
        % themselves from the beginning, now we check for the ones that
        % only map to themselves after fixing the components that have
        % a single unique cp
        % if (all_self_map)
        %     only_self_map = find((sum(Rs_Maps, 2) == 1));
        % 
        %     E = [];
        %     for j = only_self_map'
        %         B = Rs_Mats(:, :, j);
        %         C = kron(eye(MatMul), B) - kron(B', eye(MatMul));
        %         E = [stacked_fixed_sc; C]; % Append the computed matrix C to stacked_A 
        %     end
        %     [U, S, V] = svd(E);
        % 
        % 
        % end

        %% This section investigates the cyclic components
        % to what given that the unique ones must map to themselves
        Rc_Maps = zeros(3*Rc);
    
        Rc_Mats = zeros(MatMul, MatMul, 3*Rc);
        Rc_char_polys = zeros(3*Rc, MatMul+1);
        Uc = [Algs{i}{2} Algs{i}{3} Algs{i}{4}];
        for j = 1:(3*Rc)
            Rc_Mats(:, :, j) = reshape(Uc(:, j), MatMul, MatMul);
            Rc_char_polys(j, :) = charpoly(Rc_Mats(:, :, j));
        end
        
        Rc_char_polys;
    
        % Given the character polynomials of all of the columns in the
        % symmetric component, find all of the unique character polynomials
        [un_cp, ~, I] = unique(Rc_char_polys, "rows");

        % Get the count of how many times each cp shows up
        a_counts = accumarray(I,1);

        % Get the index of all of the unique cp that appear only once and store
        % it in unique_idx
        idx = find(a_counts == 1);
        Rc_unique_idx = zeros(1, length(idx));
        for j = 1:length(idx) % In case there is more than one singleton
            Rc_unique_idx(j) = find(all(Rc_char_polys == un_cp(idx(j), :), 2));
        end

        for p = 1:(3*Rc)
            if ismember(p, Rc_unique_idx)
                continue
            end
    
            same_cp_idx = find(all(Rc_char_polys == Rc_char_polys(p, :), 2));
            for q = same_cp_idx'
                Bp = Rc_Mats(:, :, p);
                Bq = Rc_Mats(:, :, q);
                D = [stacked_fixed_sc; kron(eye(MatMul), Bp) - kron(Bq', eye(MatMul))];
                
                [U, S, V] = svd(D);

                % if (p == q)
                % [p, q]
                % diag(S)
                % V(:, end-3:end)
                % end
                
                if (p == q && nnz(diag(S) > 1e-13) < (MatMul^2 - 1))
                    Rc_Maps(p, q) = 1;
                end
                if (p ~= q && nnz(diag(S) > 1e-13) < (MatMul^2))
                    Rc_Maps(p, q) = 1;
                end
                % if (nnz(diag(S) > 1e-13) < (MatMul^2 - 1))
                %     Rs_Maps(p, q) = 1;
                % end
                % Check the null space of D, if numzeros is greater than one,
                % then we store the fact that p maps to q
            end
        end
        Rc_Maps;
        has_zero_row = any(all(Rc_Maps == 0, 2));
        if (has_zero_row)
            fprintf("\t\tA cyclic component maps to nothing\n")
        else
            fprintf("\t\tAll cyclic components maps to something\n")
        end

    end
end
