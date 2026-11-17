%% Refit frequency: never vs. monthly vs. weekly vs. daily
% Same rolling-origin setup as rolling_origin_eval.m: from every day in 2014,
% forecast the next H days. The difference is that the model's parameters are
% re-estimated on all data up to the origin (an expanding window) on a fixed
% schedule:
%
% * Never   - parameters stay at their 2012-2013 estimates
% * Monthly - refit every 30 days
% * Weekly  - refit every 7 days
% * Daily   - refit every day
%
% The model structure (orders and regressors) stays the same; only the
% coefficients are updated. The normalization (training mean and std) also
% stays fixed, so all fits are on the same scale.
%
% Refitting is slow (daily = about 365 estimations per model), so only the top
% nTop models by training BIC are evaluated. Needs sarima_models.mat.

clear; clc; close all

%% Settings
H          = 28;                 % forecast up to 4 weeks ahead
schedules  = [Inf 30 7 1];       % refit every N days (Inf = never)
schedNames = ["Never" "Monthly" "Weekly" "Daily"];
nTop       = 3;                  % number of models to evaluate

%% Load
S = load('sarima_models.mat');
models = S.models;                         % {name, model, Xtrain}

yN    = [S.yTrain(:);    S.yTest(:)];      % normalized
yRaw  = [S.yTrainRaw(:); S.yTestRaw(:)];   % kW
dates = [S.timeTrain(:); S.timeTest(:)];

keep  = dates < datetime(2015,1,1);        % drop the partial final day
yN = yN(keep); yRaw = yRaw(keep); dates = dates(keep);

nTr   = numel(S.yTrain);
nAll  = numel(yN);
mu    = S.yTrainMean;
sigma = S.yTrainStd;
t0    = S.timeTrain(1);

%% Pick the top models by training BIC
nM  = size(models,1);
bic = zeros(nM,1);
for i = 1:nM
    s = summarize(models{i,2});
    bic(i) = s.BIC;
end
[~, order] = sort(bic);
useModels  = order(1:min(nTop, nM))';
nU = numel(useModels);
nS = numel(schedules);

fprintf('Evaluating %d models x %d refit schedules:\n', nU, nS);
fprintf('  %s\n', models{useModels,1});

%% Rolling forecasts with refitting
% F(t, h, model, schedule) = forecast of day t made h days earlier (kW)
F         = nan(nAll, H, nU, nS);
seconds   = zeros(nU, nS);
interPath = nan(nAll, nU);                 % intercept after each daily refit
arPath    = nan(nAll, nU);                 % AR{1} after each daily refit

for a = 1:nU
    i   = useModels(a);
    ord = ordersFromName(models{i,1});     % [p q P Q]
    X   = buildX(dates, t0, S.K, S.m, regTypeFromName(models{i,1}));

    for s = 1:nS
        mdl = models{i,2};                 % start from the 2012-2013 fit
        Fi  = nan(nAll, H);
        tic
        for o = nTr : nAll-1               % o = last observed day
            % Refit on all data up to the origin when the schedule says so.
            % (At o = nTr the data is exactly the training data, so skip.)
            if isfinite(schedules(s)) && o > nTr && mod(o - nTr, schedules(s)) == 0
                mdl = refitModel(mdl, yN(1:o), X(1:o,:), ord);
                if schedules(s) == 1
                    interPath(o,a) = mdl.Intercept;
                    if ord(1) > 0, arPath(o,a) = mdl.AR{1}; end
                end
            end

            hh = min(H, nAll - o);
            f  = forecast(mdl, hh, 'Y0', yN(1:o), 'X0', X(1:o,:), ...
                          'XF', X(o+1:o+hh,:));
            Fi(sub2ind([nAll H], o+(1:hh), 1:hh)) = f*sigma + mu;
        end
        F(:,:,a,s)   = Fi;
        seconds(a,s) = toc;
        fprintf('%-45s %-8s %6.0f s\n', models{i,1}, schedNames(s), seconds(a,s));
    end
end

% Benchmark: seasonal naive = repeat the last observed week
Fnaive = nan(nAll, H);
for o = nTr : nAll-1
    hh = min(H, nAll - o);
    lastWeek = yRaw(o-6:o);
    Fnaive(sub2ind([nAll H], o+(1:hh), 1:hh)) = lastWeek(mod(0:hh-1, 7) + 1);
end

%% Error vs. horizon for each refit schedule
common = (nTr + H) : nAll;                 % same target days for every horizon
rmse   = squeeze(sqrt(mean((yRaw(common) - F(common,:,:,:)).^2, 1, 'omitnan')));
rmse   = reshape(rmse, H, nU, nS);         % keep 3-D even if nU = 1
rmseNaive = sqrt(mean((yRaw(common) - Fnaive(common,:)).^2, 1, 'omitnan'))';

figure('Name', 'Refit schedules')
tiledlayout('flow')
for a = 1:nU
    nexttile; hold on
    plot(1:H, squeeze(rmse(:,a,:)), 'LineWidth', 1.3)
    plot(1:H, rmseNaive, 'k--', 'LineWidth', 1.5)
    xlabel('Horizon (days ahead)'); ylabel('RMSE (kW)'); grid on
    title(models{useModels(a),1}, 'Interpreter', 'none', 'FontSize', 8)
    hold off
end
legend([schedNames, "Seasonal naive"], 'Location', 'eastoutside')

% If the refit curves sit below "Never", the parameters drift over 2014 and
% refitting helps. The gain usually shows up most at long horizons, because
% short horizons are already corrected by the latest data.

%% Summary table
hShow = [1 7 14 28];  hShow = hShow(hShow <= H);
rows  = {};
for a = 1:nU
    for s = 1:nS
        r = rmse(:,a,s);
        rows(end+1,:) = {string(models{useModels(a),1}), schedNames(s), ...
                         r(hShow)', mean(r), seconds(a,s)}; %#ok<SAGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'Model','Refit','RMSE_h','MeanRMSE','Seconds'});
T = splitvars(T, 'RMSE_h', 'NewVariableNames', cellstr("h" + hShow));
T = sortrows(T, 'MeanRMSE')

fprintf('Seasonal naive mean RMSE: %.1f kW\n', mean(rmseNaive));

%% How much do the parameters drift? (daily refits)
% A steady drift in the intercept means the overall level is changing, which a
% model fit once on 2012-2013 cannot follow. A stable AR coefficient means the
% day-to-day dynamics are not changing.
if any(schedules == 1)
    figure('Name', 'Parameter drift')
    tiledlayout(2,1)
    nexttile; plot(dates, interPath, '.-'); grid on
    ylabel('Intercept (normalized)'); title('Intercept after each daily refit')
    nexttile; plot(dates, arPath, '.-'); grid on
    ylabel('AR\{1\}'); title('AR(1) coefficient after each daily refit')
    legend(string(models(useModels,1)), 'Interpreter', 'none', 'Location', 'best')
end

%% Long-horizon forecasts: never vs. daily refit (best model)
[~, iMin] = min(mean(rmse(:,:,1), 1));     % best model without refitting
testIdx   = nTr+1 : nAll;
sNever    = find(schedules == Inf, 1);
sDaily    = find(schedules == 1, 1);

if ~isempty(sNever) && ~isempty(sDaily)
    figure('Name', 'Effect of refitting'); hold on
    plot(dates(testIdx), yRaw(testIdx), 'Color', [0.6 0.6 0.6])
    plot(dates(testIdx), F(testIdx, H, iMin, sNever), 'r', 'LineWidth', 1.2)
    plot(dates(testIdx), F(testIdx, H, iMin, sDaily), 'b', 'LineWidth', 1.2)
    legend('Actual', sprintf('%d days ahead, never refit', H), ...
           sprintf('%d days ahead, refit daily', H), 'Location', 'best')
    title(models{useModels(iMin),1}, 'Interpreter', 'none')
    ylabel('kW'); hold off
end

%% Functions
function mdl = refitModel(prev, y, X, ord)
% Re-estimate a regARIMA model with the same structure, starting the optimizer
% from the previous estimates (warm start) to save time.
    p = ord(1); q = ord(2); P = ord(3); Q = ord(4);
    args = {};
    init = {'Intercept0', prev.Intercept, 'Beta0', prev.Beta, ...
            'Variance0', prev.Variance};
    if p > 0, args = [args {'ARLags',  1:p}];     init = [init {'AR0',  cell2mat(prev.AR(1:p))}];         end
    if q > 0, args = [args {'MALags',  1:q}];     init = [init {'MA0',  cell2mat(prev.MA(1:q))}];         end
    if P > 0, args = [args {'SARLags', 7*(1:P)}]; init = [init {'SAR0', cell2mat(prev.SAR(7*(1:P)))}];    end
    if Q > 0, args = [args {'SMALags', 7*(1:Q)}]; init = [init {'SMA0', cell2mat(prev.SMA(7*(1:Q)))}];    end

    tmpl = regARIMA(args{:});
    try
        mdl = estimate(tmpl, y, 'X', X, init{:}, 'Display', 'off');
    catch
        mdl = estimate(tmpl, y, 'X', X, 'Display', 'off');   % cold start
    end
end

function ord = ordersFromName(name)
% Reads [p q P Q] from a name like "... w/ (1,0,0)x(1,0,1)_7".
    tok = regexp(name, '\((\d+),0,(\d+)\)x\((\d+),0,(\d+)\)', 'tokens', 'once');
    assert(~isempty(tok), 'Could not read the orders from model name: %s', name);
    ord = str2double(tok);
end

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