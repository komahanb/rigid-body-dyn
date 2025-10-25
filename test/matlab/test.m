theta = [0.4,0.5,0.6]
C = rot(theta)

delta=1.0e-7
  omega = [1 1 1]
dt = 0.1;


dtheta1 = [0.4+delta,0.5,0.6]
dtheta2 = [0.4,0.5+delta,0.6]
dtheta3 = [0.4+delta,0.5,0.6+delta]

 CDOT =  (rot(dtheta1)-C)/delta *omega(1) + (rot(dtheta2)-C)/delta *omega(2) +  (rot(dtheta3)-C)/delta *omega(3)


