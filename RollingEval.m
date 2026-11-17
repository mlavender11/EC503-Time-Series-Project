%% Rolling-origin evaluation: how far ahead can we trust a forecast?
% From every day in 2014 (the "origin"), each model forecasts the next H days.
% Errors are grouped by horizon (1 day ahead, 2 days ahead, ...).
%
% * Error vs. horizon -> how often to reforecast.
% * Average error across horizons -> which model is best.
%
% Every horizon is scored on the SAME set of target days, so a jump in the
% curve reflects the horizon, not a change in which days are being scored.
%
% Parameters stay fixed at their training (2012-2013) estimates; only the data
% used to start each forecast is updated. Needs sarima_models.mat.

clear; clc; close all

%% Settings
H         = 28;    % forecast up to 4 weeks ahead
tol       = 0.10;  % "acceptable" = RMSE at most 10% above the 1-day RMSE
pathEvery = 14;    % in the path plots, draw a forecast path every N days

%% Load
S = load('sarima_models.mat');
models = S.models;                         % {name, model, Xtrain}
nM     = size(models,1);

yN    = [S.yTrain(:);    S.yTest(:)];      % normalized
yRaw  = [S.yTrainRaw(:); S.yTestRaw(:)];   % kW
dates = [S.timeTrain(:); S.timeTest(:)];

% The last "day" (2015-01-01) holds only a single 00:00 reading, not a full
% day of data, so drop it.
keep  = dates < datetime(2015,1,1);
yN = yN(keep); yRaw = yRaw(keep); dates = dates(keep);

nTr   = numel(S.yTrain);
nAll  = numel(yN);
mu    = S.yTrainMean;
sigma = S.yTrainStd;
t0    = S.timeTrain(1);

names = [string(models(:,1)); "Seasonal naive"];
nF    = nM + 1;

%% Rolling forecasts
% F(t, h, i) = model i's forecast of day t, made h days earlier (kW)
F = nan(nAll, H, nF);

for i = 1:nM
    X  = buildX(dates, t0, S.K, S.m, regTypeFromName(models{i,1}));
    Fi = nan(nAll, H);
    for o = nTr : nAll-1                   % o = last observed day
        hh = min(H, nAll - o);             % can't forecast past the data
        f  = forecast(models{i,2}, hh, 'Y0', yN(1:o), 'X0', X(1:o,:), ...
                      'XF', X(o+1:o+hh,:));
        Fi(sub2ind([nAll H], o+(1:hh), 1:hh)) = f*sigma + mu;
    end
    F(:,:,i) = Fi;
    fprintf('Done: %s\n', models{i,1});
end

% Benchmark: seasonal naive = repeat the last observed week
Fi = nan(nAll, H);
for o = nTr : nAll-1
    hh = min(H, nAll - o);
    lastWeek = yRaw(o-6:o);
    Fi(sub2ind([nAll H], o+(1:hh), 1:hh)) = lastWeek(mod(0:hh-1, 7) + 1);
end
F(:,:,nF) = Fi;

E = yRaw - F;                              % errors (implicit expansion)

%% Error vs. horizon (same target days for every horizon)
common = (nTr + H) : nAll;                 % days forecast at all horizons
rmse   = squeeze(sqrt(mean(E(common, :, :).^2, 1, 'omitnan')));   % H x nF

figure('Name', 'Error vs horizon'); hold on
plot(1:H, rmse(:, 1:nM), 'LineWidth', 1.2)
plot(1:H, rmse(:, end), 'k--', 'LineWidth', 2)
xlabel('Horizon (days ahead)'); ylabel('RMSE (kW)'); grid on
title('Forecast error vs. horizon (rolling origin, 2014)')
legend(names, 'Interpreter', 'none', 'Location', 'eastoutside'); hold off

%% Model comparison
hShow = [1 7 14 28];  hShow = hShow(hShow <= H);
T = array2table(rmse(hShow,:)', 'VariableNames', "h" + hShow);
T = addvars(T, names, 'Before', 1, 'NewVariableNames', 'Model');
T.MeanRMSE = mean(rmse, 1)';
T = sortrows(T, 'MeanRMSE')

%% How often to reforecast (best model)
[~, iBest] = min(mean(rmse(:,1:nM), 1));
r    = rmse(:, iBest);
hMax = find(r > (1+tol) * r(1), 1) - 1;
if isempty(hMax), hMax = H; end

fprintf('\nBest model: %s\n', names(iBest));
fprintf('RMSE: %.1f kW at 1 day, %.1f kW at %d days\n', r(1), r(end), H);
fprintf('Within %.0f%% of the 1-day error up to %d day(s) ahead\n', 100*tol, hMax);
fprintf('Beats seasonal naive at %d of %d horizons\n', sum(r < rmse(:,end)), H);

%% Forecast paths for each model
% Each colored line is one H-day forecast, starting from a different origin.
% Watch each path start near the data and relax toward the seasonal curve.
origPlot = nTr : pathEvery : nAll-1;
testIdx  = nTr+1 : nAll;

figure('Name', 'Forecast paths')
tiledlayout('flow', 'TileSpacing', 'compact')
for i = 1:nF
    nexttile; hold on
    plot(dates(testIdx), yRaw(testIdx), 'Color', [0.75 0.75 0.75])
    for o = origPlot
        hh = min(H, nAll - o);
        t  = o + (1:hh);
        plot(dates(t), F(sub2ind([nAll H nF], t, 1:hh, i*ones(1,hh))), 'LineWidth', 1.1)
    end
    title(names(i), 'Interpreter', 'none', 'FontSize', 8)
    ylabel('kW'); hold off
end

%% Forecasts at fixed horizons (best model)
% The same model's forecast of each day, made 1, 7 and 28 days in advance.
figure('Name', 'Fixed horizons'); hold on
plot(dates(testIdx), yRaw(testIdx), 'k')
plot(dates(testIdx), F(testIdx, 1, iBest),  'b')
plot(dates(testIdx), F(testIdx, 7, iBest),  'm')
plot(dates(testIdx), F(testIdx, H, iBest),  'r')
legend('Actual', '1 day ahead', '7 days ahead', sprintf('%d days ahead', H), ...
       'Location', 'best')
title("Forecasts at fixed horizons: " + names(iBest), 'Interpreter', 'none')
ylabel('kW'); hold off

%% Functions
function regType = regTypeFromName(name)
    if contains(name, 'Weekday'),     regType = "dow";
    elseif contains(name, 'Weekend'), regType = "weekend";
    elseif contains(name, 'Sunday'),  regType = "sunday";
    else,                             regType = "fourier";
    end
end

function X = buildX(dates, t0, K, m, regType)
    dates = dates(:);
    X   = fourierX(days(dates - t0), K, m);
    dow = weekday(dates);                     % 1 = Sun ... 7 = Sat
    switch regType
        case "dow"
            D = dummyvar(categorical(dow, 1:7));
            X = [X, D(:, 2:7)];               % Sunday baseline
        case "sunday"
            X = [X, double(dow == 1)];
        case "weekend"
            X = [X, double(dow == 1 | dow == 7)];
    end
end

function X = fourierX(t, K, m)
    t = t(:);
    X = zeros(numel(t), 2*K);
    for k = 1:K
        X(:, 2*k-1) = sin(2*pi*k*t/m);
        X(:, 2*k)   = cos(2*pi*k*t/m);
    end
end