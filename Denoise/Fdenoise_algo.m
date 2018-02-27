function [X_opt,h] = Fdenoise_algo(b, lambda, max_iter, lb, ub)
fprintf('***********************************\n');
fprintf('*Solving with denoise MFISTA**\n');
fprintf('***********************************\n');
fprintf('#iteration  function-value relative difference\n');
fprintf('---------------------------------------------------------------------------------------\n');

tic;
if isempty(max_iter)
    max_iter = 100;
end

[m,n] = size(b);
p = zeros(m-1,n);
q = zeros(m,n-1);

loop_tol = 1e-4;
t_pre = 1;
p_pre = p;
q_pre = q;
Gr = zeros(m,n);

for k=1:max_iter
    L_pq = L_fun(p_pre,q_pre);
    b_L = b - lambda*L_pq;
    proj_b = proj_C(b_L, ub, lb);
    [L_t_p, L_t_q] = L_tr_fun(proj_b);
    p_pr = p_pre + (1/(8*lambda))*L_t_p;
    q_pr = q_pre + (1/(8*lambda))*L_t_q;
    [p_nxt,q_nxt] = projection_pq(p_pr,q_pr);
    
    t_nxt = (1 + sqrt(1 + 4*t_pre^2))/2;
    p_pre = p_nxt + ((t_pre-1)/(t_nxt))*(p_nxt - p_pre);
    q_pre = q_nxt + ((t_pre-1)/(t_nxt))*(q_nxt - q_pre);
    t_pre = t_nxt;
    Gr_old = Gr;
    Gr = proj_b;
    h.val(k) = -norm(b_L - proj_b ,'fro')^2 + norm(b_L,'fro')^2;
    rel_dif = norm(Gr - Gr_old,'fro')/norm(Gr,'fro');
    fprintf('          %5d                 %10.10f            %10.10f \n',k,h.val(k),rel_dif);
    if k>5 && rel_dif <= loop_tol
        break;
    end
end
X_opt = proj_C(b - lambda*L_fun(p_pre,q_pre), ub, lb);
h.time = toc;
end