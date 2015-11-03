function [cdot]  = getApproxCdot(theta, theta_dot)
 % FD perturbation step size
 h=1.0d-6;  
 % approximated cdot
 cdot =  differ_rotmat(theta,h)*theta_dot(1) + ...
         differ_rotmat(theta,h)*theta_dot(2) +...
         differ_rotmat(theta,h)*theta_dot(3);
return;