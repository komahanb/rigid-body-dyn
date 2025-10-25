%% TEST differentiation of C matrix
clear; clc;

RPD = 1.0/180.0;

% these are euler angles (not the actual orientations)
theta1= 0*RPD;
theta2= 0*RPD;
theta3 = 360*RPD;

theta = [theta1 theta2 theta3];  %rad

dt    = 4;

h= 1.0d-10;

theta_dot =  theta./dt ; %[h h h];

%theta_dot = [0 0 2*RPD];            %rad/s

omega_approx = getApproxOmega(theta,theta_dot)

%which gives omega
omega = (angrate(theta)*theta_dot')'
%sdot(theta)

%(omega_approx-omega)./omega

% (1) correct the C and S matrix per book and verify with approximated
% values (YES)

% (2) See if the runs are OK
% (3) verify omega with omega_approx
% (4) Extract SDOT from -C*CDOT

%%
  
%          theta = 90*RPD;
%          C=rot([0 0 theta]);
         
%          a=[1 1 1];
%          C*a';