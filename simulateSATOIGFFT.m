function [X] = simulateSATOIGFFT(Nsim, M, flagspline, t, s, a_j, beta_j, delta_j, a, q, subCase)
%
%
% Simulation of one forward increment of S-IG subordinated BM 
% through CDF inversion via FFT. 
%
% INPUT
% M:            parameter that monitors the discretization of FFT grid
% flagspline:   1 for spline, 2 for linear interpolation (in the inversion)
% t:            yearfrac between settlement & second time instant
% s:            yearfrac between settlement & first time instant      
% a_j:          coupling parameter
% beta_j:       skew parameter 
% delta_j:      volatility
% a:            common param
% q:            self-similar parameter
% subCase:      0 in case of X subordinator, 1 in case of Z
%
% OUTPUT
% X:          simulations of values st P(X)=U (U uniform rv simulated)
%
% NB: uncomment last rows to obtain a plot of the PDF, but this takes much
% time (No FFT to compute the PDF)
%

FntNm = 'Times';
FntSz = 20;

%% Parameters
if subCase == 0
    
    beta   = 1/a_j;
    lambda = 1-a*sqrt(a_j);

elseif subCase == 1

    beta   = 1;
    lambda = a;

end

%% Assumption 2 parameters
alpha      = 0.5;
IBDaycount = 3;
B          = 1;
b          = 0.01;
omega      = 0.9;
g2         = - 1;
pplus      = -g2 - 1;  
% ashift = 0;
ashift     = 0.5*(pplus +1);

%% characteristic function of the increment of the process (between s & t)

phi_t   = @(z) exp(2*sqrt(pi)*(sqrt(beta) - sqrt(beta - 1i.*z.*t^q))*lambda);
phi_s   = @(z) exp(2*sqrt(pi)*(sqrt(beta) - sqrt(beta - 1i.*z.*s^q))*lambda);
phi_s_t = @(z) phi_t(z)./phi_s(z);
% phi_s_t = @(z) exp(2*sqrt(pi)*(- sqrt(beta - 1i.*z.*t^q) + sqrt(beta - 1i.*z.*s^q))*lambda);

% check the assumption 2
% boundF  = @(z) B.*exp(-b.*z.^omega);
% zz = [1:10000]';
% 
% figure
% plot(zz, abs(phi_s_t(zz)), 'DisplayName', 'CharFun')
% hold on
% plot(zz, boundF(zz), 'DisplayName', 'Boundary')
% grid on
% legend

%% FFT parameters
N           = 2^M;
h           = (2.*pi.*ashift./(b*N.^(omega))).^(1/(omega+1));
% h           = 0.1;
% h           = 0.7;
delta       = t-s;
D           = 8;
gammaparam  = (2*pi)/(N*h);
X0          = -gammaparam*(N-1)*0.5;
Xk          = -X0;
z1          = -0.5*h*(N-1);
zn          = -z1;

% in a struct
Params.z1   = z1;
Params.zn   = zn;
Params.dz   = h;
Params.x1   = X0;
Params.xn   = Xk;
Params.dx   = gammaparam;
Params.M    = M;
xx          = X0:gammaparam:Xk;

%% CDF on whole grid
f       = @(z) phi_s_t(z-1i.*ashift)./(1i.*z+ashift);  % f to be transformed via Fourier

%% CDF on truncated domain
point = D*sqrt(delta);
idx   = find(xx >= point);

% if last point of the grid is less than D*sqrt(delta) keep the last point
if isempty(idx)
    idx = length(xx);
end

% new grid
Xk   = xx(idx(1)).*(abs(xx(idx(1)) - point) < abs(xx(idx(1)-1) - point))+...
     xx(idx(1)-1).*(abs(xx(idx(1)) - point) > abs(xx(idx(1)-1) - point));
X0   = -Xk;
xx   = X0:gammaparam:Xk;
Phat = 1 - exp(-ashift.*xx)/(2*pi).*computeIntegral(f, xx, [], Params, 1);

% % to get rid of the obscillations
PosIdx      = find(xx>0);
FirstPosIdx = PosIdx(1);
xxNegGrid   = xx(1:FirstPosIdx);
PhatNegGrid = zeros(length(xxNegGrid), 1);
xxPosGrid   = xx(FirstPosIdx:end);
PhatPosGrid = Phat(FirstPosIdx:end);
idxPosNeg   = find(PhatPosGrid < 0);
PhatPosGrid(idxPosNeg) = 0;
PhatPosGrid(end) = 1;
xx   = [xxNegGrid, xxPosGrid]';
% 
% % put to 0 values before the CDF is increasing 
ind1       = find(PhatPosGrid>1);
PhatPosGrid(ind1) = 1;
ind2       = find(PhatPosGrid<0);
PhatPosGrid(ind2) = 0;
diff       = PhatPosGrid(2:end) - PhatPosGrid(1:end-1);
check      = find(diff<0);
% if ~isempty(check)
%     for i = 1:length(PhatPosGrid)-1
%         PhatPosGrid(i) = PhatPosGrid(i).*(all(diff(i:end)>=0)).*(PhatPosGrid(i)>=0);
%     end
%     
% end
Phat = [PhatNegGrid', PhatPosGrid]';

% plot of the CDF on the wanted strip
% idxPos  = find(xx>=0);
% xxPos   = xx(idxPos);
% PhatPos = Phat(idxPos);

% figure
% ax = gca;
% plot(xx,Phat,'LineWidth',2)
% grid on
% title('CDF via FFT on yearfraction of ', num2str(t-s), FontSize=FntSz, FontName=FntNm)
% ax.XAxis.FontSize = FntSz;
% ax.YAxis.FontSize = FntSz;
% ax.XAxis.FontName = FntNm;
% ax.YAxis.FontName = FntNm;
% xlabel('x', FontSize=FntSz, FontName=FntNm)
% ylabel('Phat',FontSize=FntSz, FontName=FntNm)

% 
% UpExcludedIdxs   = find(Phat >= 0.99999);
% UpFirstIdx       = UpExcludedIdxs(1);
% DownExcludedIdxs = find(Phat <= 1e-7);
% DownLastIdx      = DownExcludedIdxs(end);
% xx   = xx(DownLastIdx:UpExcludedIdxs);
% Phat = Phat(DownLastIdx:UpExcludedIdxs);

%% simulating & interpolating the grid values

% uncomment to simulate between max & min of the CDF (no akima needed then)
% X1 = min(PhatPosGrid);
% U  = X1 + (Phat(end-1) - X1).*rand(Nsim, 1);

U = rand(Nsim, 1); % uniform rv simulations

% interpolation
switch flagspline
    
    case 1 % spline
            [Phatunique, idxunique, ~] = unique(Phat, 'stable'); 
            xnew = xx(idxunique);
            xnew(1) = 0;
            X    = interp1(Phatunique, xnew, U, 'spline').*(interp1(Phatunique, xnew, U, 'spline')>=0);
            idx1 = find(X<xnew(1));
            idx2 = find(X>xnew(end));
            idx  = [idx1; idx2];
%             idx                        = find(X<0); % indexes where spline does fail
%             X(idx)                     = 0;

            % better to use this method when spline fails
            if ~isempty(idx)
                newX   = interp1(Phatunique, xnew, U(idx), 'makima');
                X(idx) = newX;
            end
            
    case 2 % linear
        
        [Phatunique, idxunique, ~] = unique(Phat, 'stable'); 
        X  = interp1(Phatunique, xx(idxunique), U);
end

%% PDF (uncomment to have a plot NB: a bit slow)
% density = @(x) ashift.*(exp(-ashift.*x)/pi).*h.*sum(real(exp(-1i.*([0:N-1]+0.5).*h.*x).*phi_s_t(([0:N-1]+0.5).*h-1i.*ashift)./(1i.*([0:N-1]+0.5).*h+ashift)))+...
%                         (exp(-ashift.*x)/pi).*h.*sum(real(1i*([0:N-1]+0.5).*h.*exp(-1i.*([0:N-1]+0.5).*h.*x).*phi_s_t(([0:N-1]+0.5).*h-1i.*ashift)./(1i.*([0:N-1]+0.5).*h+ashift)));
 
% adjustment as before
% densityres             = arrayfun(@(i) density(xx(i)), [1:length(xx)]');
% idx                    = find(Phat);
% densityres(1:idx(1)-1) = 0;
% 
% figure
% plot(xx, densityres, 'LineWidth', 2)
% grid on
% axis([-1.5 1.5 0 10])
% xlabel('x')
% ylabel('p')
% title('PDF on yearfraction of ', num2str(t-s))
% legend('PDF')

end