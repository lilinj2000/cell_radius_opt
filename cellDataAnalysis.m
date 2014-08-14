clc;
clear;

close all;

global LAC_CI;
global LAT_LONG;
global CELL_ANGLE;
global CELL_FREQ;
global CELL_CLASS;
global CELL_POWER;
global CUSTOM_RADIUS;
global CELL_TALIM;
global CELL_ACCMIN;
global REAL_RADIUS;

% init the cell
cell_data = initialCellData();


for ii=1:length(cell_data)
    cell_lat_long(ii, 1) = cell_data{ii, LAT_LONG}(1);
    cell_lat_long(ii, 2) = cell_data{ii, LAT_LONG}(2);
    
    real_radius(ii) = cell_data{ii, REAL_RADIUS};
end

% plot the cell geo location
figure;
plot(cell_lat_long(:, 1), cell_lat_long(:, 2), '^');
for ii = 1:length(cell_lat_long)
    if real_radius(ii)>1000
        text(cell_lat_long(ii,1), cell_lat_long(ii,2)+1, ...
            {cell_data{ii,LAC_CI}; ...
            num2str(cell_data{ii, REAL_RADIUS})});
    end
    % text(cell_lat_long(ii,1), cell_lat_long(ii,2)+1, num2str(cell_data{ii, REAL_RADIUS}));
end

figure; % radius scatter dia
% plot the real_radius
plot(real_radius, '*');
grid on;

figure; % radius bar dia
% 500, 1000, 2000, 4000, 6000, 8000, 10000, 12000
x = [500, 1000, 2000, 4000, 6000, 8000, 10000, 12000];
nbins = [250, 750, 1500, 3000, 5000, 7000, 9000, 11000];
N = hist(real_radius, nbins);
bar(x, N);
for ii = 1 : length(N)
    text(x(ii), N(ii)+5, num2str(N(ii)));
end

% % tag, freq, power, cell_type
% data_analysis = cell(length(cell_data), 4);
% 
% for ii = 1 : length(cell_data)
%     data_analysis{ii, 1} = cell_data(ii, REAL_RADIUS);
%     data_analysis{ii, 2} = cell_data{ii, CELL_FREQ};
%     data_analysis{ii, 3} = cell_data{ii, CELL_POWER};
%     data_analysis{ii, 4} = cell_data{ii, CELL_CLASS};
% end

figure;  % radius dia by cell type
hold on;
MACRO = 10;
MICRO = 20;
OTHER = 30;

for ii = 1 : length(cell_data)
    cell_type = char(cell_data{ii, CELL_CLASS});
    
    switch cell_type
        case 'MACRO'
            plot(MACRO, real_radius(ii), '*');
        case 'MICRO'
            plot(MICRO, real_radius(ii), '*');
        otherwise
            plot(OTHER, real_radius(ii), '*');
    end
end

xlim([0, OTHER+10]);
set(gca, 'xtick', [MACRO, MICRO, OTHER]);
set(gca, 'xticklabel', {'MACRO', 'MICRO', 'OTHER'});


figure; % radius dia by freq
hold on;

GSM800 = 10;
GSM900 = 20;
GSM1800 = 30;
GSM1900 = 40;
OTHER = 50;

for ii = 1 : length(cell_data)
    cell_freq = char(cell_data{ii, CELL_FREQ});
    
    switch cell_freq
        case 'GSM800'
            plot(GSM800, real_radius(ii), '*');
        case 'GSM900'
            plot(GSM900, real_radius(ii), '*');
        case 'GSM1800'
            plot(GSM1800, real_radius(ii), '*');
        case 'GSM1900'
            plot(GSM1900, real_radius(ii), '*');
        otherwise
            plot(OTHER, real_radius(ii), '*');
    end
end

xlim([0, OTHER+10]);
set(gca, 'xtick', [GSM800, GSM900, GSM1800, GSM1900, OTHER]);
set(gca, 'xticklabel', {'GSM800', 'GSM900', 'GSM1800', 'GSM1900', 'OTHER'});

figure; % radius dia by power
cell_power = [];
for ii = 1 : length(cell_data)
    cell_power(ii) = str2num(char(cell_data{ii, CELL_POWER})); 
end
plot(cell_power, real_radius, '*');
xlim([min(cell_power)-5, max(cell_power)+5]);

figure; % radius dia by freq & power
hold on;
view(3);

for ii = 1 : length(cell_data)
    cell_freq = char(cell_data{ii, CELL_FREQ});
    
    switch cell_freq
        case 'GSM800'
            plot3(GSM800, str2num(char(cell_data{ii, CELL_POWER})), real_radius(ii), '*');
        case 'GSM900'
            plot3(GSM900, str2num(char(cell_data{ii, CELL_POWER})), real_radius(ii), '*');
        case 'GSM1800'
            plot3(GSM1800, str2num(char(cell_data{ii, CELL_POWER})), real_radius(ii), '*');
        case 'GSM1900'
            plot3(GSM1900, str2num(char(cell_data{ii, CELL_POWER})), real_radius(ii), '*');
        otherwise
            plot3(OTHER, str2num(char(cell_data{ii, CELL_POWER})), real_radius(ii), '*');
    end
end

xlim([0, OTHER+10]);
set(gca, 'xtick', [GSM800, GSM900, GSM1800, GSM1900, OTHER]);
set(gca, 'xticklabel', {'GSM800', 'GSM900', 'GSM1800', 'GSM1900', 'OTHER'});