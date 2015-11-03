function [cdot]  = getApproxCdot(theta, theta_dot)
 % FD perturbation step size
 h=1.0d-10;  
 % approximated cdot
 cdot =  differ_rotmat(theta,h)*theta_dot(1) + ...
         differ_rotmat(theta,h)*theta_dot(2) +...
         differ_rotmat(theta,h)*theta_dot(3);
     
     %C3
%          C(1,1) = -sin(theta(3));
%          C(1,2) = cos(theta(3));
%          C(1,3) = 0.0d0;
%          
%          C(2,1) = -cos(theta(3));
%         C(2,2) = -sin(theta(3));
%          C(2,3) =  0.0d0;
%          
%          C(3,1) = 0.0d0;
%          C(3,2) = 0.0d0;
%          C(3,3) = 0.0d0;
         %C;
         %(cdot);
         %cdot  =C;
         %cdot_actual = C
return;