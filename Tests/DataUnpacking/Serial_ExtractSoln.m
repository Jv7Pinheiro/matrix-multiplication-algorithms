Rs4Rc1 = cell(1, 1);
Rs1Rc2 = cell(1, 1);
gd = 1;

for out = 1:2
    for ths = 1:5
        for suc = 1:78
            if CP_Data{2}{suc, ths, out}.cp_errors.rnd == 0
                Rs1Rc2{gd} = cellfun(@(x) round(x), CP_Data{2}{suc, ths, out}.out_cell, "UniformOutput",false);
                gd = gd+1;
            end
        end
    end
end