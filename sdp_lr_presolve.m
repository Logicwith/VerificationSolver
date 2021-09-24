% Solve the linear relaxation for the SDP stated in mosek_sdp.m, by dropping the semidefinite constraint

function fval = sdp_lr_presolve(H_bar, A_bar, b, A_hat, B_hat)
    
    n = size(H_bar, 1) * (size(H_bar, 1) + 1) / 2; % number of variables
    m = size(b, 1) + (size(B_hat, 1) * (size(B_hat, 1) - 1) / 2) + 1; % number of constraints
    
    % establishing model
    c = zeros(n, 1);
    idx = 1;
    for i = 1 : size(H_bar, 1)
        for j = i : size(H_bar, 1)
            if i == j
                c(idx) = H_bar(i, j);
            else
                c(idx) = H_bar(i, j) * 2;
            end
            idx = idx + 1;
        end
    end
    
    A = zeros(m, n);
    idx_1 = 1;
    for i = 1 : size(b, 1)
        idx_2 = 1;
        for j = 1 : size(A_bar{i}, 1)
            for k = j : size(A_bar{i}, 1)
                A(idx_1, idx_2) = A_bar{i}(j, k) * 2;
                idx_2 = idx_2 + 1;
            end
        end
        idx_1 = idx_1 + 1;
    end
    for i = 1 : size(A_hat, 1) - 1
        for j = i : size(A_hat, 1) - 1
            Q = A_hat(i, :)' * A_hat(j, :);
            idx_2 = 1;
            for k = 1 : size(Q, 1)
                for l = k : size(Q, 1)
                    if k == l
                        A(idx_1, idx_2) = Q(k, l);
                    else
                        A(idx_1, idx_2) = Q(k, l) + Q(l, k);
                    end
                    idx_2 = idx_2 + 1;
                end
            end
            idx_1 = idx_1 + 1;
        end
    end
    A(end, end) = 1;
    
    B = zeros(m, 1);
    B(1:size(b, 1)) = b;
    idx = size(b, 1) + 1;
    for i = 1 : size(B_hat, 1) - 1
        for j = i : size(B_hat, 1) - 1
            B(idx) = B_hat(i, j);
            idx = idx + 1;
        end
    end
    B(end) = 1;
    
    % inputing model
    model.obj = c;
    model.A = sparse(A);
    model.rhs = B;
    model.sense = '=';
    model.lb = zeros(n, 1);
    model.ub = Inf(n, 1);
    model.modelsense = 'max';
    
    % parameters
    params.OutputFlag = 0;
    
    % calling solver
    result = gurobi(model, params);
    fval = result.objval;
end

