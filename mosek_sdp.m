% Solving the dual of an SDP with the following structure:
% max  <H_bar, Y>;
% s.t. <A_bar(i), Y> = b_i,     i = 1, ..., N;
%      A_hat * Y * A_hat' = B_hat;
%      Y is positively semidefinite;    Y >= 0;    Y(end, end) = 1.
% where A_bar and A_hat has certain relation

function fval = mosek_sdp(H_bar, A_bar, b, A_hat, B_hat)
    
    % number of variables in primal prob
    n = size(H_bar, 1) * (size(H_bar, 1) + 1) / 2;
    % number of different class of constraints in dual prob
    m_1 = size(b, 1);
    m_2 = size(b, 1) * (size(b, 1) + 1) / 2;
    m_3 = n - 1;
    m_4 = 1;
    m = m_1 + m_2 + m_3 + m_4;
    
    % checking assumption that our model has to satisfy
    for i = 1 : size(A_bar, 1)
        if 2 * A_bar{i}(end, 1:(end-1)) ~= A_hat(i, 1:(end-1))
            return;
        end
    end
    
    % establishing model
    c = zeros(m, 1);
    c(1 : m_1) = -b;
    idx = m_1 + 1;
    for i = 1 : m_1
        for j = i : m_1
            c(idx) = - B_hat(i, j);
            idx = idx + 1;
        end
    end
    c(end) = 1;
    
    A = A_hat(1:end-1, 1:end-1);
    a = zeros(m_3, m);
    a(:, (m_1+m_2+1):end-1) = eye(m_3);
    idx = 1;
    for i = 1 : size(A, 2)
        for j = i : size(A, 2) + 1
            if j == size(A, 2) + 1
                a(idx , 1 : m_1) = A(:, i)';
            else
                Q = A(:, i) * A(:, j)';
                idx2 = m_1 + 1;
                for k = 1 : size(Q, 1)
                    for l = k : size(Q, 1)
                        a(idx, idx2) = Q(k, l) + Q(l, k);
                        idx2 = idx2 + 1;
                    end
                end
                if i == j
                    a(idx, (m_1+1):(m_1+m_2)) = a(idx, (m_1+1):(m_1+m_2)) / 2;
                end     
            end
            idx = idx + 1;
        end
    end
    
    idx = 1;
    blc = zeros(m_3, 1);
    for i = 1 : size(H_bar, 1) -1
        for j = i : size(H_bar, 1)
            if i == j
                blc(idx) = - H_bar(i, j);
            else
                blc(idx) = - (H_bar(i, j) + H_bar(j, i));
            end
            idx = idx + 1;
        end
    end
    buc = blc;
    
    blx = [-Inf(m_1 + m_2, 1); zeros(m_3 + m_4, 1)];
    bux = Inf(m, 1);
    
    subk = zeros(m_3, 1);
    subl = zeros(m_3, 1);
    idx = 1;
    for i = 1 : size(H_bar, 1) - 1
        for j = 1 : size(H_bar, 1) - i + 1
            subk(idx) = i + j - 1;
            subl(idx) = i;
            idx = idx + 1;
        end
    end
    
    % inputing model
    prob.c = sparse(c);
    prob.a = sparse(a);
    prob.blc = sparse(blc);
    prob.buc = sparse(buc);
    prob.blx = sparse(blx);
    prob.bux = sparse(bux);
    
    prob.bardim = size(H_bar, 1);
    prob.barc.subj = 1;
    prob.barc.subk = size(H_bar, 1);
    prob.barc.subl = size(H_bar, 1);
    prob.barc.val = 1;
    prob.bara.subi = (1:m_3)';
    prob.bara.subj = ones(m_3, 1);
    prob.bara.subk = subk;
    prob.bara.subl = subl;
    prob.bara.val = ones(m_3, 1);
    
    % parameters
    param.MSK_IPAR_NUM_THREADS = 12;
    param.MSK_IPAR_INTPNT_STARTING_POINT = 1;
    %param.MSK_IPAR_LOG = 0;
    %param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-12;
    %param.MSK_DPAR_INTPNT_CO_TOL_DFEAS = 1e-12;
    %param.MSK_DPAR_INTPNT_CO_TOL_PFEAS = 1e-12;
    
    % calling MOSEK to solve the SDP
    [~, res] = mosekopt('minimize info echo(0)', prob, param);
    fval = res.sol.itr.pobjval;
end

