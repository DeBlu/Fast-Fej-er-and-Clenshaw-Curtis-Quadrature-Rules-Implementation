function a = fejer2idst(f,N) % coefficients for Feje?r's second rule
x = cos(pi*(1:N+1)'/(N+2)); % Filippi points
fx = feval(f,x).*sin(pi*(1:N+1)'/(N+2)); % f evaluated at these points
a = idst(fx); % Chebyshev coefficients