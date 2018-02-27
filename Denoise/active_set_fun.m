function [I_not, grad_not] = active_set_fun(u, lambda, grad_u)
epl = 10^-5;

I = (u == -lambda & grad_u > epl) |...
    (u == lambda & grad_u < -epl);
I_not = not(I);
grad_not = grad_u(I_not);
end