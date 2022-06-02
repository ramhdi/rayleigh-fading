%% Jake's method for Rayleigh fading channel model generation
%% Adapted from: http://www.ee.iitm.ac.in/~giri/pdfs/EE5141/Jakes-Simulation.pdf 
%% ramhdi, 26/11/2020

function r = fading2(len, fd, Ts)
  N = 100; %34, the more the better
  a = 2*pi*rand(1,N); % random phase, uniform distrib. 0 to 2*pi
  b = 2*pi*rand(1,N); % random phase, uniform distrib. 0 to 2*pi
  alpha = 2*pi*rand(1,N); % random angle of incidence, uniform distrib. 0 to 2*pi
  
  t = (1:len)*Ts; % samples
  
  ri = zeros(1,len); rq = zeros(1,len);
  for m = 1:N
    ri = ri + 1/sqrt(N) * cos(2*pi*fd*cos(alpha(m))*t + a(m));
    rq = rq + 1/sqrt(N) * cos(2*pi*fd*cos(alpha(m))*t + b(m));
  end
  
  r = ri + j*rq;