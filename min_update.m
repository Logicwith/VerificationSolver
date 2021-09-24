% Solve the minimum for our quadratic program and refresh the status of our problem

function [status, min] = min_update(H, p, A, b, c)
    %[~, min, exf] = cplexqp(2 * H, 2 * p, [], [], A, b, zeros(1, size(A, 2)));
    
    [~, min, exf] = gurobiqp(H, p, A, b);
    
    if exf == -2
        status = 3;
        return;
    end
    if c < min
        status = 2;
        return;
    end
    status = 0;
end