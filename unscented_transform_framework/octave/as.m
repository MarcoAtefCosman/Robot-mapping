clc,clear all,close all;
w_c=[1,2,3,4,5];
sigma_points=[1,2,3,4,5;1,2,3,4,5];
mu=[55;55];
n=length(mu);
sigma = zeros(n,n);

for i=1:2*n+1
    sigma = sigma + w_c(1,i)*bsxfun(@minus,sigma_points,mu)*bsxfun(@minus,sigma_points,mu)';
endfor