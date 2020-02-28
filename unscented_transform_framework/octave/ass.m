n = 2;
sigma_points=[1,2,3,4,5;1,2,3,4,5];
mu=[55;55];
w_c=[1,2,3,4,5];
sigma = zeros(n,n);
for i=1:2*n+1
    vare=bsxfun(@minus,sigma_points,mu)*bsxfun(@minus,sigma_points,mu)'
    sigma = sigma + w_c(1,i)*vare
    w_c(1,i)
endfor