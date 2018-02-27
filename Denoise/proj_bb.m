function x = proj_bb(u,lambda)
x = min(max(u,-ones(size(u))*lambda),ones(size(u))*lambda);
end