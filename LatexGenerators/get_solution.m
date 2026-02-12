%% Rs1 Rc16

for i = 1:6
    fprintf("\\subsection{Rs1\\_Rc16: Solution %d}\n", i)

    factor_matrices = cyc2fac(Rs1_Rc16_Solutions{i});
    string = "Rs1-Rc16-Solution-" + i;
    
    print_char_poly(factor_matrices, 1, string);
    fprintf("\n\\newpage\n\n")
    fprintf("\\footnotesize{\n")
    make_matrices_latex(factor_matrices, 49, 1, string);
    fprintf("}\n\n")

    fprintf("\n\\begin{landscape}\n")
    fprintf("\\resizebox{\\columnwidth}{!}{\n")

    make_incidence_latex(Rs1_Rc16_Solutions{i})
    
    fprintf("}\n")
    fprintf("\n\\vfill\n")

    fprintf("\\resizebox{\\columnwidth}{!}{\n")
    make_pairing_latex(Rs1_Rc16_Solutions{i}, 1, 16);
    
    fprintf("}\n")
    fprintf("\\end{landscape}\n")
end

%% Rs1 Rc16 - Strassen

for i = 1:1
    fprintf("\\subsection{Rs1\\_Rc16: Strassen}\n")

    factor_matrices = cyc2fac(Rs1_Rc16_Strassen{i});
    string = "Rs1-Rc16-Strassen-" + i;
    
    print_char_poly(factor_matrices, 1, string);
    fprintf("\n\\newpage\n\n")
    fprintf("\\footnotesize{\n")
    make_matrices_latex(factor_matrices, 49, 1, string);
    fprintf("}\n\n")

    fprintf("\n\\begin{landscape}\n")
    fprintf("\\resizebox{\\columnwidth}{!}{\n")

    make_incidence_latex(Rs1_Rc16_Strassen{i})
    
    fprintf("}\n")
    fprintf("\n\\vfill\n")

    fprintf("\\resizebox{\\columnwidth}{!}{\n")
    make_pairing_latex(Rs1_Rc16_Strassen{i}, 1, 16);
    
    fprintf("}\n")
    fprintf("\\end{landscape}\n")
end

%% Rs4 Rc 15 - Strassen
for i = 1:1
    fprintf("\\subsection{Rs4\\_Rc15: Solution %d}\n", i)

    factor_matrices = cyc2fac(Rs4_Rc15_Solutions{i});
    string = "Rs4-Rc15-Solution-" + i;
    
    print_char_poly(factor_matrices, 4, string);
    fprintf("\n\\newpage\n\n")
    fprintf("\\footnotesize{\n")
    make_matrices_latex(factor_matrices, 49, 4, string);
    fprintf("}\n\n")

    fprintf("\n\\begin{landscape}\n")
    fprintf("\\resizebox{\\columnwidth}{!}{\n")

    make_incidence_latex(Rs4_Rc15_Solutions{i})
    
    fprintf("}\n")
    fprintf("\n\\vfill\n")

    fprintf("\\resizebox{\\columnwidth}{!}{\n")
    make_pairing_latex(Rs4_Rc15_Solutions{i}, 4, 15);
    
    fprintf("}\n")
    fprintf("\\end{landscape}\n")
end

%% Rs13 Rc 12 - Transformed
for i = 1:1
    fprintf("\\subsection{Rs12\\_Rc12 (Transformed): Solution %d}\n", i)

    factor_matrices = cyc2fac(Rs12_Rc12_Solution_Transformed{i});
    string = "Rs12-Rc12-Solution-" + i;
    
    print_char_poly(factor_matrices, 13, string);
    fprintf("\n\\newpage\n\n")
    fprintf("\\footnotesize{\n")
    make_matrices_latex(factor_matrices, 49, 13, string);
    fprintf("}\n\n")

    fprintf("\n\\begin{landscape}\n")
    fprintf("\\resizebox{\\columnwidth}{!}{\n")

    make_incidence_latex(Rs12_Rc12_Solution_Transformed{i})
    
    fprintf("}\n")
    fprintf("\n\\vfill\n")

    fprintf("\\resizebox{\\columnwidth}{!}{\n")
    make_pairing_latex(Rs12_Rc12_Solution_Transformed{i}, 13, 12);
    
    fprintf("}\n")
    fprintf("\\end{landscape}\n")
end

%% Rs13 Rc 12
for i = 1:1
    fprintf("\\subsection{Rs13\\_Rc12: Solution %d}\n", i)

    factor_matrices = cyc2fac(Rs13_Rc12_Solutions{i});
    string = "Rs13-Rc12-Solution-" + i;
    
    print_char_poly(factor_matrices, 13, string);
    fprintf("\n\\newpage\n\n")
    fprintf("\\footnotesize{\n")
    make_matrices_latex(factor_matrices, 49, 13, string);
    fprintf("}\n\n")

    fprintf("\n\\begin{landscape}\n")
    fprintf("\\resizebox{\\columnwidth}{!}{\n")

    make_incidence_latex(Rs13_Rc12_Solutions{i})
    
    fprintf("}\n")
    fprintf("\n\\vfill\n")

    fprintf("\\resizebox{\\columnwidth}{!}{\n")
    make_pairing_latex(Rs13_Rc12_Solutions{i}, 13, 12);
    
    fprintf("}\n")
    fprintf("\\end{landscape}\n")
end


%% Rs16 Rc11

for i = 1:32
    fprintf("\\subsection{Rs16\\_Rc11: Solution %d}\n", i)

    factor_matrices = cyc2fac(Rs16_Rc11_Solutions{i});
    string = "Rs16-Rc11-Solution-" + i;
    
    print_char_poly(factor_matrices, 16, string);
    fprintf("\n\\newpage\n\n")
    fprintf("\\footnotesize{\n")
    make_matrices_latex(factor_matrices, 49, 16, string);
    fprintf("}\n\n")

    fprintf("\n\\begin{landscape}\n")
    fprintf("\\resizebox{\\columnwidth}{!}{\n")

    make_incidence_latex(Rs16_Rc11_Solutions{i})
    
    fprintf("}\n")
    fprintf("\n\\vfill\n")

    fprintf("\\resizebox{\\columnwidth}{!}{\n")
    make_pairing_latex(Rs16_Rc11_Solutions{i}, 16, 11);
    
    fprintf("}\n")
    fprintf("\\end{landscape}\n")
end

% Rs16 Rc 11 - Strassen
for i = 1:1
    fprintf("\\subsection{Rs16\\_Rc11: Actual Strassen %d}\n", i)

    factor_matrices = cyc2fac(Soln);
    string = "Rs16-Rc11-Strassen-" + i;
    
    print_char_poly(factor_matrices, 16, string);
    fprintf("\n\\newpage\n\n")
    fprintf("\\footnotesize{\n")
    make_matrices_latex(factor_matrices, 49, 16, string);
    fprintf("}\n\n")

    fprintf("\n\\begin{landscape}\n")
    fprintf("\\resizebox{\\columnwidth}{!}{\n")

    make_incidence_latex(Soln)
    
    fprintf("}\n")
    fprintf("\n\\vfill\n")

    fprintf("\\resizebox{\\columnwidth}{!}{\n")
    make_pairing_latex(Soln, 16, 11);
    
    fprintf("}\n")
    fprintf("\\end{landscape}\n")
end