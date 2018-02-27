function u = TV_BB(D, y, lambda, tol)

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
    dual_gap = @(x, u)(lambda*(norm(D*x,1)) - u'*D*x);
    phi_fun = @(u)(0.5*(norm((D'*u),2))^2 - u'*grad_out);
%% Initial solve
    % alternate solution
     
    tic;
    R = chol(Hess,'lower');
    optsL.LT = true;
    u_inter = linsolve(R, grad_out, optsL);
    optsU.UT = true;
    u.chol = linsolve(R', u_inter, optsU);
    u.time_u = toc;
    
    
    % Faster solution
    %tic;
    %u.opt = D'\y; 
    %u.time2 = toc;
    %u_dual = u.opt;
    u_dual = u.chol;
    x.opt = y - D'*u_dual;
    x_dual = x.opt;
    
    if norm(u_dual, 'inf') <= lambda
        u.x = x_dual;
        u.time = toc;
        return;
    end
    
%% Full solve
    u.init = proj_u(u_dual, lambda);
    u_solve = u.init;
    t = 1;
    dl_gap = 1000;
    while dl_gap > tol
        gr_u = grad_u(u_solve);
        [I_not, grad_not] = active_set_fun(u_solve, lambda, gr_u);
        dir_grad = zeros(size(gr_u));
        dir_grad(I_not) = grad_not;
        if mod(t,2)==1
            al = norm(dir_grad,2)^2/dot(dir_grad,Hess*dir_grad);
        else
            al = dot(dir_grad,Hess*dir_grad)/norm(Hess*dir_grad,2)^2;
        end
        u_solve = proj_bb(u_solve - al*dir_grad, lambda);
        x_dual = y - D'*u_solve;
        t = t + 1;
        dl_gap = dual_gap(x_dual, u_solve);
        fprintf('%.5d Duality gap \n', dl_gap);
    end
    u.x = x_dual;
    u.time = toc;
end