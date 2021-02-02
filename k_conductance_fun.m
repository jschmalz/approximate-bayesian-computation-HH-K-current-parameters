function [t,g] = k_conductance_fun(V,theta_p)

% parameters
theta_n = theta_p(1:end-1);
g_k = theta_p(end);


% initial conditions and run parameters
N = size(V,2);
n0 = (0.24/g_k)^(1/4);
tspan = [0:0.1:12];
nt = size(tspan,2);
n = zeros(nt,N);

% run for each constant membrane potential
for i = 1:N
    [t,ni] = ode45(@(t,n) odeK(t,n,theta_n,V(i)), tspan, n0);
    n(:,i) = ni; 
end

% calculate conductance
g = g_k*n.^4;

end