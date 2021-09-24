% Verifies whether a value c can be attained
% by a convex quadratic function f(x) = x'*H*x + 2*p'*x 
% on a polyhedral {x : Ax = b, x >= 0}.
% We assume that H = U*D*U' for some U, D with positive entries,
% and A, b, p has positive entries.

% eps: parameter controlling the cutting plane.
% Possible values of status:
% 0 - not sure;
% 1 - attined;
% 2 - not attained;
% 3- infeasible.

% Requires MOSEK toolbox to run.

function [status, minval, maxval_ub, maxval_lb] = verification_solver(H, D, U, p, A, b, c, eps)

    % calculating minimum and status
    
    [status, minval] = min_update(H, p, A, b, c);
    if status ~= 0
        maxval_ub = -Inf;
        maxval_lb = -Inf;
        return;
    end
    
    % initializing
    
    lb_global = -Inf;
    ub_global = Inf;
    
    % performing vertex searching, SDR update and cutting plane methods
    
    while check_feasibility(A, b) ~= 0       
        [x_large, lb_1] = search_large_vertex(H, U, p, A, b);
        [x_lmax, lb_2, updated] = search_local_max_vertex(H, p, A, b, x_large);
        lb_current = max([lb_1; lb_2]);
        ub_current = max_ub_sdr_update(H, p, A, b, c);
        ub_global = max([lb_global; min([ub_current; ub_global])]);
        lb_global = max([lb_current; lb_global]);
        
        if c > ub_global || c < lb_global
            status = 2;
            break;
        end
        
        if updated == 1
            cut1 = konno_cut(H, p, A, b, x_lmax, lb_global, eps);
            cut2 = konno_cut(H, p, A, b, x_large, lb_global, eps);
            %cut3 = konno_cut(y, lb, eps);
            cuts = [cut1 cut2];
        else
            cut1 = konno_cut(H, p, A, b, x_large, lb_global, eps);
            %cut3 = konno_cut(y, lb, eps);
            cuts = cut1;
        end
        
        [H, U, p, A, b] = update_feasible_region(U, D, p, A, b, cuts);
    end
    
    % returning variables
    
    maxval_ub = ub_global;
    maxval_lb = lb_global;
end

