%% Numerically differentiate the rotation matrix using chain rule
function [pC_ptheta] = differ_rotmat(theta, dtheta)
%%PDC
pC_ptheta = (rot(theta+dtheta) - rot(theta) )/dtheta;