len = length(Rs16_Rc11_Solutions);


for i = 1:len
    KX = Rs16_Rc11_Solutions{i};
    FKX = cyc2fac(KX);
    for j = 1:len
        if i ~= j
            KY = Rs16_Rc11_Solutions{j};
            FKY = cyc2fac(KY);

            bool = EquivalentUpToPermutation(FKX{1}, FKX{2}, FKX{3}, FKY{1}, FKY{2}, FKY{3});
            if bool == 1
                fprintf('Solution %d and %d result in true\n', i, j);
            end

        end
    end
end

fprintf('Done\n');



% soln1 = 5;
% soln2 = 11;
% 
% 
% KX = Rs16_Rc11_Solutions{soln1};
% KY = Rs16_Rc11_Solutions{soln2};
% 
% FKX = cyc2fac(KX);
% FKY = cyc2fac(KY);
% 
% equiv_up_to_perm(FKX{1}, FKX{2}, FKX{3}, FKY{1}, FKY{2}, FKY{3})