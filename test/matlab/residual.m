function [res] = residual(m, g0, q, qdot,t)
                  
r           = q(1:3);
theta       = q(4:6);
v           = q(7:9);
omega       = q(10:12);

r_dot       = qdot(1:3);
theta_dot   = qdot(4:6);
v_dot       = qdot(7:9);
omega_dot   = qdot(10:12);

c = r*m;

J = -skew(r)*skew(r)*m;
      
C = rot(theta);
      
S = angrate(theta);

res(1:3)    =    C*r_dot' - v';
res(4:6)    =    S*theta_dot' - omega';
res(7:9)    =    m*v_dot - cross(c,omega_dot) + cross(omega, (m*v - cross(c,omega))) -  sin(2*t); %-m*g0;
res(10:12)  =    cross(c, v_dot)' + J*omega_dot' + skew(c)*skew(omega)*v' + skew(omega)*J*omega' -J*omega_dot';

end