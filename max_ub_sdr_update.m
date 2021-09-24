% Update the upper bound for the maximum value via SDR.

function [ub, x] = max_ub_sdr_update(H, p, A, b, c)

    % setting parameters in SDP
    H_bar = [H p; p' 0];
    A_bar = cell(size(A, 1), 1);
    for i = 1 : size(A, 1)
        A_bar{i} = [zeros(size(A, 2)) 0.5*A(i, :)'; 0.5*A(i, :) 0];
    end
    A_hat = [A zeros(size(A, 1), 1); zeros(1, size(A, 2)) 0];
    B = b * b';
    B_hat = [B zeros(size(B, 1), 1); zeros(1, size(B, 2)) 0];
    %A_h = [A'*A zeros(size(A, 2), 1); zeros(1, size(A, 2)) 0];
    
    % presolve the linear relaxation of SDP
    if exist('c', 'var')
        ub = sdp_lr_presolve(H_bar, A_bar, b, A_hat, B_hat);
        if ub < c
            return;
        end
    end
    
    % calling SDP solver
    ub = mosek_sdp(H_bar, A_bar, b, A_hat, B_hat);
    
end