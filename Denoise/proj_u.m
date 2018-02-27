function u_pr = proj_u(u, lambda)

u_pr = (-lambda<=u & u<=lambda).*u + (u>lambda).*lambda + (u<-lambda).*-lambda;

end