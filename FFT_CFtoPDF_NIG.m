function PdfNIG = FFT_CFtoPDF_NIG(params)

% Inversion of the Characteristic function via FFT with NIG

%% Characteristic function

% model parameters
delta   = params(1); % 
gamma   = params(2); % 
beta    = params(3); % 
CharFun = @(u) exp(- delta * (sqrt(gamma^2 - (beta + 1i.*u).^2) - (sqrt(gamma^2 - beta^2))));

%% Discretization grid
Npow = 17;
N    = 2^Npow; % discretization points
A    = 2000;   % domain upper boundary

% integral domain grid
eta  = A/N;               % spacing on the domain grid
v    = [0:eta:A*(N-1)/N]; % int domain grid
v(1) = 1e-22;             % adjust starting point near zero
CharFunVals = CharFun(v);

% density grid
lambda = 2 * pi/(N * eta);                  % relation
k      = - lambda * N/2 + lambda * (0:N-1); % density grid

% 
w   = ones(1, N); w(1) = 0.5; w(end) = 0.5;                % trapezoidal rule weights
l   = w .* eta .* CharFunVals .* exp( 1i * pi * (0:N-1) ); % function in the DFT
PdfNIG = real( fft(l) / pi );                              % apply discrete Fourier transform


end