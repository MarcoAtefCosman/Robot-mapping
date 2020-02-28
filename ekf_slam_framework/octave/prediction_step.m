function [mu, sigma] = prediction_step(mu, sigma, u)
% Updates the belief concerning the robot pose according to the motion model,
% mu: 2N+3 x 1 vector representing the state mean
% sigma: 2N+3 x 2N+3 covariance matrix
% u: odometry reading (r1, t, r2)
% Use u.r1, u.t, and u.r2 to access the rotation and translation values

% TODO: Compute the new mu based on the noise-free (odometry-based) motion model
% Remember to normalize theta after the update (hint: use the function normalize_angle available in tools)

%the new mu is also has the size of 2N+3*1 where only the first 3 elements are effected.
%we need to define a control vector u and a transform matrix F
distance=u.t;
rotation1=u.r1;
rotation2=u.r2;
old_theta=mu(3)

U=[distance*cos(rotation1+old_theta);distance*sin(rotation1+old_theta);normalize_angle(rotation1+rotation2)];
F=[eye(3);zeros(length(mu)-3,3)];

mu=mu+F*U;
mu(3) = normalize_angle(mu(3));
% TODO: Compute the 3x3 Jacobian Gx of the motion model
%Gxt=[f1_x,f1_y,f1_theta;f2_x,f2_y,f2_theta;f3_x,f3_y,f3_theta]
Gxt=[0,0,-distance*sin(old_theta+rotation1);0,0,distance*cos(old_theta+rotation1);0,0,0];

% TODO: Construct the full Jacobian G
G=eye(size(sigma))+F*Gxt*F';

% Motion noise
motionNoise = 0.1;
R3 = [motionNoise, 0, 0; 
     0, motionNoise, 0; 
     0, 0, motionNoise/10];
R = zeros(size(sigma,1));
R(1:3,1:3) = R3;

% TODO: Compute the predicted sigma after incorporating the motion
sigma=G*sigma*G'+R;

end
