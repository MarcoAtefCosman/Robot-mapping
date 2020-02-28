function [mu, sigma, sigma_points] = prediction_step(mu, sigma, u)
% Updates the belief concerning the robot pose according to the motion model.
% mu: state vector containing robot pose and poses of landmarks obeserved so far
% Current robot pose = mu(1:3)
% Note that the landmark poses in mu are stacked in the order by which they were observed
% sigma: the covariance matrix of the system.
% u: odometry reading (r1, t, r2)
% Use u.r1, u.t, and u.r2 to access the rotation and translation values

% For computing lambda.
global scale;

% Compute sigma points
sigma_points = compute_sigma_points(mu, sigma);

% Dimensionality
n = length(mu);
% lambda
lambda = scale - n;

% TODO: Transform all sigma points according to the odometry command
% Remember to vectorize your operations and normalize angles
% Tip: the function normalize_angle also works on a vector (row) of angle
for i=1:2*n+1
  sigma_points(1,i)=sigma_points(1,i)+u.t*cos(sigma_points(3,i)+u.r1);
  sigma_points(2,i)=sigma_points(2,i)+u.t*sin(sigma_points(3,i)+u.r1);
  sigma_points(3,i)=sigma_points(3,i)+(u.r1+u.r2);
  sigma_points(3,i)=normalize_angle(sigma_points(3,i));
endfor

% Computing the weights for recovering the mean
wm = [lambda/scale, repmat(1/(2*scale),1,2*n)];
wc = wm;

% TODO: recover mu.
% Be careful when computing the robot's orientation (sum up the sines and
% cosines and recover the 'average' angle via atan2)
x=0;
y=0;
mu=zeros(n,1);

% TODO: recover mu.
% Be careful when computing the robot's orientation (sum up the sines and
% cosines and recover the 'average' angle via atan2)
for i=1:2*n+1
    mu=mu+wm(i)*sigma_points(:,i);
    x=x+wm(i)*cos(sigma_points(3,i));
    y=y+wm(i)*sin(sigma_points(3,i));    
endfor
mu(3)=atan2(y,x);
mu(3)=normalize_angle(mu(3));

% TODO: Recover sigma. Again, normalize the angular difference
sigma=zeros(n);
for i=1:2*n+1
    difference=sigma_points(:,i)-mu;
    difference(3)=normalize_angle(difference(3));
    sigma=sigma+wc(i)*difference*difference';
endfor

% Motion noise
motionNoise = 0.1;
R3 = [motionNoise, 0, 0; 
     0, motionNoise, 0; 
     0, 0, motionNoise/10];
R = zeros(size(sigma,1));
R(1:3,1:3) = R3;

% TODO: Add motion noise to sigma
sigma=sigma+R;

end