%% Load Data
client = 'MT_210';
[xh,y] = load_electricity(client);

%% Weekly Naive
snWeekly = lag(xh(:, client), 168).(client);
testPeriod = xh.Time >= datetime(2014, 1, 1) & xh.Time < datetime(2015, 1, 1);
validSamples = testPeriod & ~isnan(snWeekly) & ~isnan(y);

yTest = y(validSamples);
yTestPred = snWeekly(validSamples);
tT = xh.Time(validSamples);
e = yTest - yTestPred;

ME = mean(e);
MAE = mean(abs(e));
RMSE = sqrt(mean(e.^2));
nMAE = MAE / mean(yTest) * 100;

%%

wkStart = datetime(2014,3,3);                  % try a few different weeks
w = tT >= wkStart & tT < wkStart + days(7);

figure
plot(tT(w), yTest(w), 'k', tT(w), yTestPred(w), 'r--')
legend('Actual', 'Seasonal naive (weekly)'); ylabel('kWh'); grid on


E = timetable(tT, e, abs(e), 'VariableNames', {'Err','AbsErr'});
Ed = retime(E, 'daily', 'mean');

figure
plot(Ed.tT, Ed.AbsErr); ylabel('Daily MAE (kWh)'); grid on

byHour = groupsummary(E, 'tT', 'hourofday', 'mean', {'Err','AbsErr'});
byDay  = groupsummary(E, 'tT', 'dayname',   'mean', {'Err','AbsErr'});
byMonth  = groupsummary(E, 'tT', 'monthname', 'mean', {'Err','AbsErr'});

figure
histogram(e, 100); xlabel('e (kWh)'); ylabel('Count')

results = table(ME, MAE, RMSE, nMAE, ...
    'RowNames', {'SeasonalNaiveWeekly'})