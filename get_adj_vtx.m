% Compute an adjacent vertex from a given original vertex
% A, b: description of the polytope
% basic_pos: the position of basic variables of the original vertex
% enter_idx: the (enter_idx -th)'s variable in the nonbasic variable list will become the entering variable

function adj_vtx = get_adj_vtx(A, b, basic_pos, enter_idx)
    non_basic_pos = sort(setdiff(1:size(A, 2), basic_pos));
    pivot_col_pos = non_basic_pos(enter_idx);
    
    % minimal ratio test
    B = A(:, basic_pos);
    f = B \ b;
    d = B \ A(:, pivot_col_pos);
    min_ratio = inf;
    for i = 1 : size(A, 1)
        if d(i) > 0
            ratio = f(i) / d(i);
            if ratio < min_ratio
                min_ratio = ratio;
                pivot_row_pos = i;
            end
        end
    end
    
    % entering and leaving the basis
    if pivot_row_pos ~= 0
        new_basic_pos = basic_pos;
        new_basic_pos(pivot_row_pos) = pivot_col_pos;
        new_basic_pos = sort(new_basic_pos);
        B = A(:, new_basic_pos);
        adj_vtx = zeros(size(A, 2), 1);
        adj_vtx(new_basic_pos) = B \ b;
    end
end

