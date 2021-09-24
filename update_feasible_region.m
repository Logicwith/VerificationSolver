% Update the feasible region by the given cutting planes

function [H, U, p, A, b] = update_feasible_region(U_0, D_0, p_0, A_0, b_0, cuts)
    n = size(cuts, 2);
    A = [A_0  zeros(size(A_0, 1), n); cuts'  -eye(n)];
    b = [b_0; ones(n, 1)];
    p = [p_0; zeros(n, 1)];
    U = [U_0; zeros(n, size(U_0, 2))];
    H = U * D_0 * U';
end

