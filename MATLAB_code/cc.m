function a = cc(f,N) % (N+1)-coefficients for C-C quadrature
x = cos(pi*(0:N)'/N); % C-C points
fx = feval(f,x)/(2*N); % f evaluated at these points
g = fft(fx([1:N+1 N:-1:2])); % FFT
a = [g(1); g(2:N)+g(2*N:-1:N+2); g(N+1)]; % Chebyshev coefficients