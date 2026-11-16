client = 'MT_210';
[xh,y] = load_electricity(client);

%% Weekly Plots
nWeeks = 10;
firstDay = dateshift(xh.Time(1), 'start', 'day', 'next'); % Get first monday
allDays = firstDay:days(1):(xh.Time(end) - days(7)); % Array of days from first monday to last monday
mondays = allDays(weekday(allDays) == 2); % Extract all mondays

pick = sort(mondays(randperm(numel(mondays), nWeeks))); % Select 10 random mondays

figure; hold on
for k = 1:nWeeks
    wk = xh(timerange(pick(k), pick(k) + days(7)), :);
    hoursIn = hours(wk.Time - pick(k));     % 0 to 167 hours since Monday 00:00
    plot(hoursIn, wk.(client))
end
hold off
xticks(0:24:168)
xticklabels({'Mon','Tue','Wed','Thu','Fri','Sat','Sun',''})
xlabel('Day of week'); ylabel('kWh per hour'); grid on
title(sprintf('Weekly Plots of %s for Randomly Chosen Weeks', client));
ylabel('Energy (kWh per hour)')
