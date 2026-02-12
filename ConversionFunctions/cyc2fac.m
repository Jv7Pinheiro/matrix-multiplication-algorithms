function T = cyc2fac(D)
%CYC2FAC Converts four cylic matrices into three factor matrices
%   
%   T = cyc2fac({S U V W}), returns the ktensor T:
%       T{1} = {S U V W}
%       T{2} = {S W U V}
%       T{3} = {S V W U}
%
    A = [D{1} D{2} D{3} D{4}];
    B = [D{1} D{4} D{2} D{3}];
    C = [D{1} D{3} D{4} D{2}];
    T = {A B C};
end