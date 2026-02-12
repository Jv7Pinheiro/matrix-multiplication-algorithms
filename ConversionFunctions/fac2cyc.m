function [S, U, V, W] = fac2cyc(T, Rs, Rc)
    S = T{1}(:, 1:Rs);
    Ru = Rs+1
    Rv = Rc+Ru-1
    U = T{1}(:, Ru:Rv);
    Rv = Rv+1
    Rw = Rv+Rc-1
    V = T{1}(:, Rv:Rw);
    Rw = Rw+1
    W = T{1}(:, Rw:end);
end