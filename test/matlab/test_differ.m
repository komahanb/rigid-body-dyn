%% TEST differentiation of C matrix
clear; clc;

RPD = 1.0/180.0;

theta1= 20*RPD;
theta2= 15*RPD;
theta3= 10*RPD;

theta = [theta1 theta2 theta3];  %rad
theta_dot = [1 1 .1];            %rad/s

omega_approx = getApproxOmega(theta,theta_dot)

%which gives omega
omega = (angrate(theta)*theta_dot')'
%sdot(theta)

%(omega_approx-omega)./omega