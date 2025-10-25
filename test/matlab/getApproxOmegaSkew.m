function [omega_skew] = getApproxOmegaSkew(theta, theta_dot)
 
cdot_approx =getApproxCdot(theta,theta_dot);

c = rot(theta);

omega_skew = -cdot_approx*c';