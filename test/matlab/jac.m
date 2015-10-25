function [jac] =jac(a, m, q, qdot)

r           = q(1:3);
theta       = q(4:6);
v           = q(7:9);
omega       = q(10:12);

r_dot       = qdot(1:3);
theta_dot   = qdot(4:6);
%v_dot       = qdot(7:9);
%omega_dot   = qdot(10:12);

c = r*m;

J = -skew(r)*skew(r)*m;
      
C = rot(theta);
      
S = angrate(theta);

S_dot = sdot(theta);

jac(1:3,1:3)        = a*C;
jac(1:3,4:6)        = skew(C*r_dot')*S;
jac(1:3,7:9)        = -eye(3,3);
jac(1:3,10:12)      = zeros(3,3);

jac(4:6,1:3)        = zeros(3,3);
jac(4:6,4:6)        = S_dot + skew(S*theta_dot')*S + a*S;
jac(4:6,7:9)        = zeros(3,3);
jac(4:6,10:12)      = -eye(3,3);

jac(7:9,1:3)        = zeros(3,3);
jac(7:9,4:6)        = zeros(3,3);
jac(7:9,7:9)        = m*(a*eye(3,3) + skew(omega));
jac(7:9,10:12)      = -a*skew(c) + skew(skew(c)*omega')-m*skew(v) -skew(omega)*skew(c);


jac(10:12,1:3)      = zeros(3,3);
jac(10:12,4:6)      = zeros(3,3);
jac(10:12,7:9)      = a*skew(c) + skew(c)*skew(omega);
jac(10:12,10:12)    = a*J - skew(c)*skew(v) +  skew(omega)*J -skew(J*omega');

end