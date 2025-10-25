% octave script to check the correctness of linear system setup
  m=5.0d0 
   r       = [ 1,2,3 ]
   theta   = [ 10.0*pi/180, (-20.0d0)*pi/180, (-40.0d0)*pi/180 ]
   v       = [ 1.0d0, -2.0d0, 3.0d0 ]
   omega   = [ -1.0d0, 3.0d0, -4.0d0 ]

   Y       = [r, theta, v, omega ]

   r_dot       = [ 1.0d0, -2.0d0, 3.0d0 ]
   theta_dot   = [ -1.0d0, 2.0d0, -1.0d0 ]
   v_dot       = [ 1.0d0, 2.0d0, -2.0d0 ]
   omega_dot   = [ -1.0d0, 4.0d0, -2.0d0 ]

   YPRIME       = [r_dot, theta_dot, v_dot, omega_dot ]


   g0 = [ -1.0d0, -1.0d0, -1.0d0]
   re = [ -1.0d0, 2.0d0, 1.0d0]

   c = re*m

   J = [25, 10,5;
          10 10 -10;
          5 -10 25];

C=[0.757858  -0.633226  0.157083;
    0.558673   0.754219  0.345019;
    -0.336950  -0.173717  0.925358]

  S          =  [ 0.76586  -0.63323  0.00000;
                   0.64300   0.75422  0.00000;
                   0.00000  -0.17372  1.00000 ]

  C*r_dot' - v'
                  
  S*theta_dot' - omega'

  m*v_dot - cross(c,omega_dot) + cross(omega, (m*v - cross(c,omega))) -m*g0

  cross(c, v_dot)' + J*omega_dot' + skew(c)*skew(omega)*v' + skew(omega)*J*omega' -J*omega_dot'

  
%% verify JACOBIAN

S_dot      =   [0.00000  0.00000  0.00000
                0.00000  0.00000  0.00000
                0.00000  0.00000  0.00000]

a=2.0 

a*C

skew(C*r_dot')*S

S_dot + skew(S*theta_dot')*S + a*S

m*(a*ones(3,3) + skew(omega))

-a*skew(c) + skew(skew(c)*omega')-m*skew(v) -skew(omega)*skew(c)

a*skew(c) + skew(c)*skew(omega)

  a*J - skew(c)*skew(v) +  skew(omega)*J -skew(J*omega')
  function [s] = skew(a)
  s =[0,-a(3) , a(2);
  a(3), 0, -a(1);
  -a(2), a(1),0];          
  end
