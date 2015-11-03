function [res] = residual(m, g0, re, q, qdot, t)
                  
r           = q(1:3);   % position of body axis from origin
theta       = q(4:6);   % orientation of body axis with inertial
v           = q(7:9);   % velocity of the body axis
omega       = q(10:12); % angular velocity of the body axis

r_dot       = qdot(1:3); %v;
theta_dot   = qdot(4:6); %omega;
v_dot       = qdot(7:9);
omega_dot   = qdot(10:12);

c = re*m;
J = -skew(re)*skew(re)*m;
      
C = rot(theta); 
%S = angrate(theta);

res(1:3)    =    C*r_dot' - v';
%res(4:6)    =    S*theta_dot' - omega';
res(4:6)    =    getApproxOmega(theta,theta_dot)- omega';
res(7:9)    =    m*v_dot - cross(c,omega_dot) + cross(omega, (m*v - cross(c,omega)));
res(10:12)  =    cross(c, v_dot)' + J*omega_dot' + skew(c)*skew(omega)*v' + skew(omega)*J*omega';

end