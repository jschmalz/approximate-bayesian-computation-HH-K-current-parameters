function dndt = odeK(t,n,theta_p,V)
% reseting membrane potential of 0

% parameters
a1 = theta_p(1); 
a2 = theta_p(2);
a3 = theta_p(3);
b1 = theta_p(4);
b2 = theta_p(5);

% alpha beta functions
alpha_n = a1*(V+a2)/(exp( (V+a2)/a3 ) -1);
beta_n = b1*exp(V/b2);

% ode
dndt = alpha_n * (1-n) - beta_n*n;
end