function a = fejer1(f,N)
x = cos(pi*(2*(0:N)'+1)/(2*N+2)); % coefficients for Feje?r's first rule
fx = feval(f,x)/(N+1); % first kind of Chebyshev points
g = fft(fx([1:N+1 N+1:-1:1])); % f evaluated at these points
hx = real(exp(2*1i*pi*(0:2*N+1)/(4*N+4)).*g'); % FFT
a = hx(1:N+1);a(1)=0.5*a(1); % Chebyshev coefficients