% Generate a Konno's cut at a given vertex x
% lb: current best fval found on vertices
% eps: parameter

function cut = konno_cut(H, p, A, b, x, bestobjval, eps)
    % get the position of basic variables for x_0
    basic_pos = zeros(size(A, 1), 1);
    idx = 1;
    for i = 1 : size(x, 1)
        if x(i) ~= 0
            basic_pos(idx) = i;
            idx = idx + 1;
        end
    end

    % degenerate case
    if basic_pos(end) == 0
        return;
    end
    
    % positions of nonbasic variables
    non_basic_pos = setdiff(1:size(A, 2), basic_pos);
    non_basic_pos = sort(non_basic_pos);
    
    % computing parameters
    B = A(:, basic_pos);
    N = A(:, non_basic_pos);
    F = B \ N;
    f = B \ b;
    H_BB = H(basic_pos, basic_pos);
    H_BN = H(basic_pos, non_basic_pos);
    H_NN = H(non_basic_pos, non_basic_pos);
    p_B = p(basic_pos);
    p_N = p(non_basic_pos);
    D = H_NN + F' * H_BB * F - 2 * F' * H_BN;
    d = p_N - F' * p_B - F' * H_BB * f + H_BN' * f;
    phi_0 = f' * H_BB * f + 2 * p_B' * f;
    l = size(A, 2) - size(A, 1);
    
    % generate Tuy's cut
    theta = zeros(l, 1);
    for i = 1 : l
        theta(i) = max(roots([D(i, i) 2*d(i) phi_0-bestobjval-eps]));
    end
    
    % generate Konno's cut
    tau = zeros(l, 1);
    for i = 1 : l
        u = [-d; bestobjval - phi_0 + eps];
        t = ones(l, 1) ./ theta;
        A_ineq = [F -f; -t' 1];
        b_ineq = zeros(size(A_ineq, 1), 1);
        A_eq = [D(i, :) d(i)];
        b_eq = 1;
        %[~, fval, exf] = cplexlp(u, A_ineq, b_ineq, A_eq, b_eq, zeros(l + 1, 1));
        [~, fval, exf] = gurobilp('min', u, A_eq, b_eq, A_ineq, b_ineq);
        if exf == 1
            tau(i) = fval;
        else
            tau(i) = theta(i);
        end
    end
    
    cut = zeros(size(x, 1), 1);
    cut(non_basic_pos) = 1 ./ tau;
end