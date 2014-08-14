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

% % plot the cell geo location
% figure;
% plot(cell_lat_long(:, 1), cell_lat_long(:, 2), '^');

figure; % plot the cell geo location
hold on;

radius_section = [0, 500, 1000, 2000, 4000, 6000, 8000, 10000];
colors = {'b^', 'g^', 'r^', 'b*', 'g*', 'r*', 'b+', 'g+'};
legend_string = {'500', '1000', '2000', '4000', '6000', '8000', ...
    '10000', '> 10000'};
legend_h = zeros(1,length(legend_string));

% plot(cell_lat_long(:, 1), cell_lat_long(:, 2), '^');
for ii = 1:length(cell_lat_long)
    
    for jj=length(radius_section) : -1 : 1
        if real_radius(ii)> radius_section(jj)
            h = plot(cell_lat_long(ii, 1), cell_lat_long(ii, 2), colors{jj});
            
            if legend_h(jj)==0
                legend_h(jj) = h;
            end
            break;
        end
    end
end

legend(legend_h, legend_string);
title('cell geo-location');

figure; % radius scatter dia
% plot the real_radius
plot(real_radius, '*');
grid on;
title('radius scatter');

figure; % radius bar dia
% 500, 1000, 2000, 4000, 6000, 8000, 10000, 12000
x_lable = {'500', '1000', '2000', '4000', '6000', '8000', '10000', '>10000'};
nbins = [250, 750, 1500, 3000, 5000, 7000, 9000, 11000];
N = hist(real_radius, nbins);
bar([1 : length(nbins)], N);
for ii = 1 : length(N)
    text(ii, N(ii)+5, num2str(N(ii)));
end
set(gca, 'xtick', [1 : length(nbins)]);
set(gca, 'xticklabel', x_lable);
title('radius distribution');

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
title('radius vs cell type');

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
title('radius vs frequence');

figure; % radius dia by power
cell_power = [];
for ii = 1 : length(cell_data)
    cell_power(ii) = str2num(char(cell_data{ii, CELL_POWER})); 
end
plot(cell_power, real_radius, '*');
xlim([min(cell_power)-5, max(cell_power)+5]);
title('radius vs power');

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
title('radius vs frequence vs power');