function [p,q] = L_tr_fun(X)

[m,n] = size(X);
p = zeros(m-1,n);
q = zeros(m,n-1);

for i=1:m-1
    for j=1:n-1
        p(i,j) = X(i,j) - X(i+1,j);
        q(i,j) = X(i,j) - X(i,j+1);
    end
end

for i=1:m-1
   p(i,end) =  X(i,end) - X(i+1,end);
end

for j=1:n-1
   q(end,j) =  X(end,j) - X(end,j+1);
end

end