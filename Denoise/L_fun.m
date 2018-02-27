function L_pq = L_fun(p,q)

m = size(q,1);
n = size(p,2);
L_pq = zeros(m,n);

p_aug = [zeros(1,n);p;zeros(1,n)];
p_aug = [zeros(size(p_aug,1),1) p_aug];

q_aug = [zeros(m,1) q zeros(m,1)];
q_aug = [zeros(1,size(q_aug,2));q_aug];

for i=2:m+1
    for j=2:n+1
        L_pq(i-1,j-1) = p_aug(i,j) + q_aug(i,j) - p_aug(i-1,j) - q_aug(i,j-1);
    end
end

end