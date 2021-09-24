% Search for a locally maximum vertex, starting from a given vertex x_0

function [x_lmax, fval, updated] = search_local_max_vertex(H, p, A, b, x_0)
    x = x_0;
    updated = 0; % whether x_lmax differs from x_0
    
    while true
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
        
        non_basic_pos = sort(setdiff(1:size(A, 2), basic_pos));
        
        % find a KKT vertex
        fx = x' * H * x + 2 * p' * x;
        is_kkt_pt = 1;
        u = 2 * H * x + 2 * p;
        u_B = u(basic_pos);
        u_N = u(non_basic_pos);
        B = A(:, basic_pos);
        N = A(:, non_basic_pos);
        r = ((B' \ u_B)' * N - u_N')';
        for i = 1 : size(non_basic_pos, 2)
            if r(i) < 0
                y = get_adj_vtx(A, b, basic_pos, i);
                is_kkt_pt = 0;
                x = y;
                updated = 1;
                break;
            end
        end
        
        % check whether the KKT vertex is locally maximum
        if is_kkt_pt == 1
            is_local_max = 1;
            for i = 1 : size(non_basic_pos, 2)
                y = get_adj_vtx(A, b, basic_pos, i);
                fy = y' * H * y + 2 * p' * y;
                if fy > fx
                    is_local_max = 0;
                    x = y;
                    updated = 1;
                    break;
                end
            end
            if is_local_max == 1
                fval = x' * H * x + 2 * p' * x;
                break;
            end
        end
    end
    
    x_lmax = x;
end