function [mu, sigma] = recover_gaussian(sigma_points, w_m, w_c)
% This function computes the recovered Gaussian distribution (mu and sigma)
% given the sigma points (size: nx2n+1) and their weights w_m and w_c:
% w_m = [w_m_0, ..., w_m_2n], w_c = [w_c_0, ..., w_c_2n].
% The weight vectors are each 1x2n+1 in size,
% where n is the dimensionality of the distribution.

% Try to vectorize your operations as much as possible

% TODO: compute mu
%mu is n,2n+1 matrix should return n,1:
mu=(w_m*sigma_points')';

% TODO: compute sigma
w_c=[w_c;w_c];
sigma=zeros(length(mu));
for i=1:length(w_c)
  sigma=sigma+w_c(i).*(sigma_points-mu)*(sigma_points-mu)';
endfor
