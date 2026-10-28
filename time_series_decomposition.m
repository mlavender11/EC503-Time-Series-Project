clear; clc; close all;
rng(1);

%% Generate example series
N = 240; % 20 yrs of monthly data
t = (1:N)';
s = 12; % Yearly seasonal period

% Generate trend, cycle, seasonal, and noise components
genTrend = 50 + 0.15*t + 0.0003*t.^2; % quadratic growth over 20 yrs
genCycle = 4*sin(2*pi*t/96) .* (1 + 0.3*sin(2*pi*t/200)); % Varying ~8 yr cycle
genSeason = 6*sin(2*pi*t/s) + 2*cos(4*pi*t/s); % Yearly cycle
genNoise = 1.5*randn(N, 1); % random white noise

y = genTrend + genCycle + genSeason + genNoise; % generated signal


%%  Estimate comps
% [LT, estSeason, estNoise] = trenddecomp(y, 'stl', s);
% [p, ~, mu]  = polyfit(t, LT, 2);
% estTrend  = polyval(p, t, [], mu); % Estimated trend
% estCycle  = LT - estTrend;

[~, ST] = trenddecomp(y, "stl", s);
yDS = y - ST;                       % contains trend + cycle + noise

w   = [0.5; ones(s-1,1); 0.5] / s; 
TC  = conv(yDS, w, 'same');
TC([1:s/2, end-s/2+1:end]) = NaN;

ok = ~isnan(TC);
[estTrend, estCycle] = hpfilter(TC(ok), Smoothing=129600);   % Econometrics Toolbox

%% Plot

figure(Name='Time Series Decomposition', Position=[100 100 800 900]);
tl = tiledlayout(5, 1, TileSpacing="compact");
title(tl, 'Time Series Decomposition')

nexttile; plot(t, y, 'k'); hold on; plot(t, estTrend+estCycle, 'r', LineWidth=1.5);
title('Original Series and Trend-Cycle'); legend('y', 'T+C', Location='best');

nexttile; plot(t, genTrend, 'k--'); hold on; plot(t, estTrend, 'r');
title('Estimated vs Generated Trend'); legend('generated', 'estimated', Location='best');

nexttile; plot(t, genCycle, 'k--'); hold on; plot(t, estCycle, 'r');
title('Estimated vs Generated Cycle'); legend('generated', 'estimated', Location='best');

nexttile; plot(t, genSeason, 'k--'); hold on; plot(t, estSeason, 'r');
title('Estimated vs Generated Seasonality'); legend('generated', 'estimated', Location='best');

nexttile; stem(t, estNoise);
title('Residual noise'); xlabel('Month')