%%
function [res] = residual(m, g0, q, qdot)
                  
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
      
%S = angrate(theta);

res(1:3)    =    C*r_dot' - v';
%res(4:6)    =    S*theta_dot' - omega';
res(4:6)    =    getApproxOmega(theta,theta_dot)- omega';
res(7:9)    =    m*v_dot - cross(c,omega_dot) + cross(omega, (m*v - cross(c,omega))) -sin(2*t);
res(10:12)  =    cross(c, v_dot)' + J*omega_dot' + skew(c)*skew(omega)*v' + skew(omega)*J*omega' -J*omega_dot';

end 
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
function [s] = skew(a)
 s =[0,    -a(3), a(2);
     a(3),  0,  -a(1);
    -a(2), a(1), 0];          
end
function [CBI] = rot(theta)

    CBI(1,1) =  cos(theta(2))*cos(theta(3)) + sin(theta(1))*sin(theta(2))*sin(theta(3));
    CBI(2,1) = -cos(theta(2))*sin(theta(3)) + sin(theta(1))*sin(theta(2))*cos(theta(3));
    CBI(3,1) =  cos(theta(1))*sin(theta(2));

    CBI(1,2) =  cos(theta(1))*sin(theta(3));
    CBI(2,2) =  cos(theta(1))*cos(theta(3));
    CBI(3,2) = -sin(theta(1));

    CBI(1,3) = -sin(theta(2))*cos(theta(3)) + sin(theta(1))*cos(theta(2))*sin(theta(3));
    CBI(2,3) =  sin(theta(2))*sin(theta(3)) + sin(theta(1))*cos(theta(2))*cos(theta(3));
    CBI(3,3) =  cos(theta(1))*cos(theta(2));
    
end
function [SIB] = angrate(theta)

    SIB(1,1) = cos(theta(3));
    SIB(2,1) = -sin(theta(3));
    SIB(3,1) = 0.0;

    SIB(1,2) = cos(theta(1))*sin(theta(3));
    SIB(2,2) = cos(theta(1))*cos(theta(3));
    SIB(3,2) = -sin(theta(1)) ;

    SIB(1,3) = 0.0;
    SIB(2,3) = 0.0;
    SIB(3,3) = 1.0;
    
end
function [SDOT] = sdot(theta)

    SDOT(1,1) = -sin(theta(3));
    SDOT(2,1) = -cos(theta(3));
    SDOT(3,1) = 0.0;

    SDOT(1,2) = -sin(theta(1))*sin(theta(3)) + cos(theta(1))*cos(theta(3));
    SDOT(2,2) = -sin(theta(1))*cos(theta(3)) - cos(theta(1))*sin(theta(3));
    SDOT(3,2) = -cos(theta(1)) ;

    SDOT(1,3) = 0.0;
    SDOT(2,3) = 0.0;
    SDOT(3,3) = 0.0;
    
end