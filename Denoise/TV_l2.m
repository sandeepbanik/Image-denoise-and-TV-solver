function u = TV_l2(D, y, lambda, tol)
%% Parameter check
    tic;
    Hess = D*D';
    grad_out = D*y;
    n = size(y, 1);
    if size(D, 2) ~= n || size(D, 1) ~= n-1
        error('The differencing matrix size is incorrect');
    end
    if isempty(tol)
        tol=1e-5;
    end
    grad_u = @(u)(Hess*u - grad_out);
    dual_gap = @(x, u)(lambda*(norm(D*x,2)) - u'*D*x);
%% Initial solve

    tic;
    R = chol(Hess,'lower');
    optsL.LT = true;
    u_inter = linsolve(R, grad_out, optsL);
    optsU.UT = true;
    u.chol = linsolve(R', u_inter, optsU);
    u.time_u = toc;

    u_solve = u.chol;
    x.opt = y - D'*u_solve;
    x_dual = x.opt;
    u_norm = norm(u_solve, 2);
    if u_norm <= lambda
        u.x = x_dual;
        u.time = toc;
        return;
    end

%% Full solve
    al_GP = 1/4;
    al_MSN = 0;
    d_gp = 1000;
    d_msn = 1000;
    if lambda < norm(y,2)
        while abs(d_gp) > tol
            u_solve = proj_l2(u_solve - al_GP*grad_u(u_solve), lambda);
            x_dual = y - D'*u_solve;
            d_gp = dual_gap(x_dual, u_solve);
            fprintf('%d dual gap \n', d_gp);
        end
    else
        while (d_msn > tol || u_norm >= lambda + 1e-4)
            R = chol(D*D' + al_MSN*eye(n-1), 'lower');
            u_temp = linsolve(R, grad_out, optsL);
            u_msn = linsolve(R', u_temp, optsU);
            q = linsolve(R', u_msn, optsU);
            u_norm = norm(u_msn,2);
            al_MSN = al_MSN - (u_norm/norm(q,2))*(1 - u_norm/lambda);
            x_dual = y - D'*u_msn;
            d_msn = dual_gap(x_dual, u_msn);
            fprintf('%d dual gap %d norm \n', d_msn, u_norm);
        end
    end
    u.x = x_dual;
    u.time = toc;
end