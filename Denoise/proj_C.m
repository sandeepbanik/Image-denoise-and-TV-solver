function proj_x = proj_C(X, ub, lb)

[m,n] = size(X);
proj_x = zeros(size(X));
for i=1:m
    for j=1:n
        if (X(i,j) <= ub) && (lb <= X(i,j))
            proj_x(i,j) = X(i,j);
        elseif (X(i,j) < lb)
            proj_x(i,j) = lb;
        elseif (ub < X(i,j))
            proj_x(i,j) = ub;
        end
    end
end

end