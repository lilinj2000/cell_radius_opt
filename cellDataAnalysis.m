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


rogers = 1;
turkey = 2;
india = 3;

rogers_area_rural_orangeville = 1;
rogers_area_suburban_brampton = 2;
rogers_area_suburban_scarborough = 3;
rogers_area_suburban = 4;
rogers_area_urban_dt = 5;
rogers_area_all = 6;

subhead = {'(rogers rural orangeville)'};

% init the cell
cell_data = initialCellData(rogers, rogers_area_rural_orangeville);


for ii=1:length(cell_data)
    cell_lat_long(ii, 1) = cell_data{ii, LAT_LONG}(1);
    cell_lat_long(ii, 2) = cell_data{ii, LAT_LONG}(2);
    
    real_radius(ii) = cell_data{ii, REAL_RADIUS};
end

radiusAnalysis(real_radius, cell_lat_long, subhead)

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
title(vertcat({'radius vs cell type'}, subhead));

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
title(vertcat({'radius vs frequence'}, subhead));

figure; % radius dia by power
cell_power = [];
for ii = 1 : length(cell_data)
    cell_power(ii) = str2num(char(cell_data{ii, CELL_POWER})); 
end
plot(cell_power, real_radius, '*');
xlim([min(cell_power)-5, max(cell_power)+5]);
title(vertcat({'radius vs power'}, subhead));

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
grid on;
title(vertcat({'radius vs frequence vs power'}, subhead));