function [omega] = getApproxOmega(theta, theta_dot)
 
cdot_approx =getApproxCdot(theta,theta_dot);

c = rot(theta);

omega_skew = -cdot_approx*c';

omega = unskew(omega_skew);