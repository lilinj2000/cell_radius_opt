% clc;
% clear;
% 
close all;

importSysVar;



% init the cell
% cell_data = initialCellData(china, china_wuxi);
% cell_data = turkey_cell;
% cell_radius = t_real_radius;

% rogers
% subhead = {'(urban dt - real radius)'};
% cell_data = urban_dt_cell;
% cell_radius = urban_dt_radius;
subhead = {'(rogers all - real radius)'};
cell_data = rogers_cell;
cell_radius = t_real_radius;


idx_real_radius = 3;  % 0.95 measurements

index_radius = 0;
cell_data_tag = [];
for ii=1:size(cell_data, 1)
    if cell_data{ii, REAL_RADIUS}~=0
        cell_data_tag = [cell_data_tag; cell_data(ii, :)];
        
        index_radius = index_radius + 1;
        cell_lat_long(index_radius, 1) = cell_data{ii, LAT_LONG}(1);
        cell_lat_long(index_radius, 2) = cell_data{ii, LAT_LONG}(2);
    
        real_radius(index_radius) = cell_radius(ii, idx_real_radius);
    end
end

radiusAnalysis(real_radius, cell_lat_long, subhead)

figure;  % radius dia by cell type
hold on;
MACRO = 10;
MICRO = 20;
OTHER = 30;

for ii = 1 : length(cell_data_tag)
    cell_type = char(cell_data_tag{ii, CELL_CLASS});
    
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

for ii = 1 : length(cell_data_tag)
    cell_freq = char(cell_data_tag{ii, CELL_FREQ});
    
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
for ii = 1 : length(cell_data_tag)
    if isempty(cell_data{ii, CELL_POWER})
        cell_power(ii) = 0;
    else
        cell_power(ii) = str2num(char(cell_data{ii, CELL_POWER})); 
        
        if ~isempty(cell_data{ii, CELL_ACCMIN})
            cell_power(ii) = cell_power(ii) + str2num(char(cell_data{ii, CELL_ACCMIN}));
        end
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

    for ii = 1 : length(cell_data_tag)
        cell_freq = char(cell_data_tag{ii, CELL_FREQ});

        switch cell_freq
            case 'GSM800'
                plot3(GSM800, cell_power(ii), real_radius(ii), '*');
            case 'GSM900'
                plot3(GSM900, cell_power(ii), real_radius(ii), '*');
            case 'GSM1800'
                plot3(GSM1800, cell_power(ii), real_radius(ii), '*');
            case 'GSM1900'
                plot3(GSM1900, cell_power(ii), real_radius(ii), '*');
            otherwise
                plot3(OTHER, cell_power(ii), real_radius(ii), '*');
        end
    end

    xlim([0, OTHER+10]);
    set(gca, 'xtick', [GSM800, GSM900, GSM1800, GSM1900, OTHER]);
    set(gca, 'xticklabel', {'GSM800', 'GSM900', 'GSM1800', 'GSM1900', 'OTHER'});
    grid on;
    title(vertcat({'radius vs frequence vs power'}, subhead));
end