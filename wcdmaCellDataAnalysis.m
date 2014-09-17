clc;
clear;

close all;

china = 1;

rogers_area_rural_orangeville = 1;
rogers_area_suburban_brampton = 2;
rogers_area_suburban_scarborough = 3;
rogers_area_suburban = 4;
rogers_area_urban_dt = 5;
rogers_area_all = 6;

turkey_beyoglu = 7;
turkey_incek = 8;
turkey_sariyer = 9;
turkey_ulus = 10;
turkey_uskudar = 11;
turkey_all = 12;

india_medium = 13;
india_rural = 14;
india_rural_highway = 15;
india_all = 16;

% china_beijing = 17;

subhead = {'(rogers urban_dt - real radius)'};

% init the cell
cell_data = initialCellData(rogers, rogers_area_urban_dt);


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


cell_power = [];
for ii = 1 : length(cell_data)
    if isempty(cell_data{ii, CELL_POWER})
        cell_power(ii) = 0;
    else
        cell_power(ii) = str2num(char(cell_data{ii, CELL_POWER})); 
    end
end

if ~isempty(find(cell_power~=0))
    figure; % radius dia by power
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
end