%%
clear;clc;close all;

% initial guess and other parameters
m           = 10.0;
q(1:12)     = 1.0;
qdot(1:12)  = 1.0;
g0          = [1,-1,1];

% setup the time marching
itmax        = 300;     % max iteration to drive residual to zero
update_tol   = 1.0e-10; % norm(dq) 
tmax         = 5.00;    % final time
dt           = 0.1;      % time step 

a           = 1.0/dt;
state(12,1000) = 0;
cnt    = 0;
time   = 0;
while (time<=tmax) 
    cnt=cnt+1;
    time = time + dt; 
    % drive the residual to zero
    for i = 1:itmax
        % get the rhs
        r = residual(m,g0,q,qdot,time)';
        % get the lhs
        J = jac(a, m, q, qdot);
        % get the update
        [dq, FLAG, ITER] = lsqr(J,-r);dq=dq';
        % check if the update is small enough
        if (norm(dq) < update_tol) 
            % stop the update(steady state reached)
            break;
        else
            % update the state
            q = q + dq;
            qdot = qdot + a *dq;   
        end
    end
    
    % save the states and time for plotting
    state(cnt,1:12) = q;
    tt(cnt) = time;
end

plotstates(state,tt);