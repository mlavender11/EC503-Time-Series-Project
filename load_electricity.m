function [xh,y] = load_electricity(client)
% Load Data
opts = detectImportOptions('LD2011_2014.txt', ...
    'Delimiter', ';', 'DecimalSeparator', ',');
opts = setvartype(opts, 1, 'datetime');
opts = setvaropts(opts, 1, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
T = readtable('LD2011_2014.txt', opts);
T.Properties.VariableNames{1} = 'Time';

TT = table2timetable(T, 'RowTimes', 'Time');

% Choose and load client

x = TT(:, client);

% Trim leading zeros (before client connected)
first = find(x.(client) > 0, 1); 
x = x(first:end, :);

% Combine 15 min readings into hourly readings
xh = retime(x, 'hourly', @(v) sum(v)/4);   % kWh per hour
y = xh.(client);
end