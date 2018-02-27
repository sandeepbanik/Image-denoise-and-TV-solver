function [r, s] = projection_pq(p, q)

m = size(q,1);
n = size(p,2);

r = zeros(size(p));
s = zeros(size(q));

for i=1:m-1
    for j=1:n-1
        r(i,j) = p(i,j)/max(1,norm([p(i,j);q(i,j)],2));
        s(i,j) = q(i,j)/max(1,norm([p(i,j);q(i,j)],2));
    end
end

for i=1:m-1
   r(i,end) =  p(i,end)/max(1,abs(p(i,end)));
end

for j=1:n-1
   s(end,j) =  q(end,j)/max(1,abs(q(end,j)));
end

end