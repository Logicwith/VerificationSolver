% Feasibility test for the region {Ax = b, x >= 0}

function [is_feasible] = check_feasibility(A, b)
    is_feasible = 0;
    [~, ~, exf] = gurobilp('min', zeros(size(A, 2), 1), A, b, [], []);
    if exf == 1
        is_feasible = 1;
    end
end

