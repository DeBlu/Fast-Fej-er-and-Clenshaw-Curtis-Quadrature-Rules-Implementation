function a = fejer1dct(f,N) % coefficients for Feje?r's second rule
x = cos(pi*(2*(0:N)'+1)/(2*N+2)); % The first kind of Chebyshev points
fx = feval(f,x); % f evaluated at these points
a = dct(fx)*sqrt(2/(N+1));a(1)=a(1)/sqrt(2); % Chebyshev coefficients