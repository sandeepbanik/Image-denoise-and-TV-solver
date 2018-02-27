function u = proj_l2(u_solve, lambda)
% x_norm = norm(u_solve,2);
% u = (x_norm>lambda)*u_solve./x_norm + u_solve*(x_norm <= lambda);
% u = (x_norm>=lambda)*(1 - lambda/x_norm).*u_solve + zeros(size(u_solve))*(x_norm < lambda);
u = min((lambda / norm(u_solve, 2)), 1) * u_solve;
end