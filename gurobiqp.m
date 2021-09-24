% Call Gurobi to solve the following quadratic program:
% max  x'*H*x + 2*p'*x
% s.t. Ax = b, x >= 0

function [x, fval, exf] = gurobiqp(H, p, A, b)
    model.Q = sparse(H);
    model.obj = 2 * p;
    model.A = sparse(A);
    model.rhs = b;
    model.sense = '=';
    model.lb = zeros(size(H, 1), 1);
    model.ub = Inf(size(H, 1), 1);
    
    params.OutputFlag = 0;
    
    result = gurobi(model, params);
    x = result.x;
    fval = result.objval;
    if strcmp(result.status, 'OPTIMAL') %#ok<*BDSCA>
        exf = 1;
    end
    if strcmp(result.status, 'INFEASIBLE') || strcmp(result.status, 'INF_OR_UNBD') || strcmp(result.status, 'UNBOUNDED')
        exf = -2;
    end
end

