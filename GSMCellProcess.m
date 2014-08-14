
close all;
clear;
clc;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  Parameter Confgure %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% General parametr configure

% Specifies the distance between one site and another when searching 
% for co-located cells. 
% If the distance is less than the value of Adjustment_Clearance, 
% the sites are to be converged into one site.
% Unit: meters
Adjustment_Clearance = 10; 

% Specifies a cell type to use when Cell Radius Adjuster fails 
% to get the cell type from insufficient information.
% Value range: Lower or Upper
Lower = 1;
Upper = 2;
CellType_Default = Lower; 

% Specifies the default value of cell radius 
% when no cell radius is calculated.
% Unit: meters
Radius_DefaultCellRadius = 500; % unit meters

% Specifies a radius for fixed cells, 
% such as PICO cell in GSM and femto cell in WCDMA.
% Unit: meters
Radius_FixedCell = 100;

% Cell type parameter configure

% Specifies the method Cell Radius Adjuster uses to get the cell type.
% Value range: PropagationModel or Frequency
PropagationModel = 1;
Frequency = 2;
CellType_JudgingMethod = Frequency;

% Specifies a value of cell boundary in propagation model. 
% If the value is greater than 1000, the cell type is judged as upper cell. 
% If the value is less than or equal to 1000, the cell type is judged 
% as lower cell.
CellType_BoundaryPropagationModel = 1000;

% Specifies a list of frequencies of upper cell.
% Example:
% UPPER_FREQUENCE_LIST =  100, 200, 10663, 10661
NA = [];
UPPER_FREQUENCE_LIST = NA;

% Specifies a list of frequencies of lower cell.
% Example: LOWER_FREQUENCE_LIST =  10, 20, 30, 50
LOWER_FREQUENCE_LIST = NA;

% ISD & Voronoi Algorithm configure

% Specifies the algorithm Cell Radius Adjuster uses to determine the cell
% radius.
% value range: ISD or Voronoi
ISD = 1;
Voronoi = 2;
Radius_AdjustmentAlgorithm = ISD;

% Specifies the scope of radius adjustment. 
% If the scope is set to Fully, the cell radius of all cells is adjusted. 
% If the scope is set to Partly, only those cells of which cell radius, 
% eutranCellPosition, or eutranCellPolygon does not exist are adjusted. 
% value range: Fully or Partly
Fully = 1;
Partly = 2;
Radius_AdjustmentBehavior = Fully;

% Specifies the boundary cell radius when Voronoi diagram is used.
Voronoi_BoundaryCellRadius = 1500;

% Specifies the factor that is configured to make cell coverage larger 
% than that is calculated from the Voronoi diagram.
Voronoi_RadiusFactor = 1.2;

% Specifies the search scope when ISD algorithm is used.
% Unit: meters
ISD_SearchScope = 4500;

% Specifies the number of closest cells chosen to calculate the radius 
% when ISD algorithm is used. 
% Example:
% ISD_NumberOfClosestCells = 5
ISD_NumberOfClosestCells = 5;

% Specifies the distance ratio that is used to adjust the distance 
% towards radius when ISD algorithm is used.
ISD_GeographicalDistanceRatio = 0.7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%Algorithm Selection%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
false = 0;
true = 1;

custom_radius = false;
ta_radius = true;
oh_propagation_model = true;
isd_radius = true;
voronoi_radius_max = true;
voronoi_radius_mean = true;

methods = 0;
index_custom_radius = 0;
index_ta_radius = 0;
index_oh_propagation_model = 0;
index_isd_radius = 0;
index_voronoi_radius_max = 0;
index_voronoi_radius_mean = 0;

methods_name = {};

if custom_radius
    methods = methods + 1;
    index_custom_radius = methods;
    methods_name{index_custom_radius} = 'custom radius';
end

if ta_radius 
    methods = methods + 1;
    index_ta_radius = methods;
    methods_name{index_ta_radius} = 'TA';
end

if oh_propagation_model
    methods = methods + 1;
    index_oh_propagation_model = methods;
    methods_name{index_oh_propagation_model} = 'OH';
end

if isd_radius
    methods = methods + 1;
    index_isd_radius = methods;
    methods_name{index_isd_radius} = 'ISD';
end

if voronoi_radius_max
    methods = methods + 1;
    index_voronoi_radius_max = methods;
    methods_name{index_voronoi_radius_max} = 'VORONOI max';
end

if voronoi_radius_mean
    methods = methods + 1;
    index_voronoi_radius_mean = methods;
    methods_name{index_voronoi_radius_mean} = 'VORONOI mean';
end

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

cell_data = initialCellData();

% [ m, 1]
radius_inner = zeros(length(cell_data), 1);

% the real radius is [ m, 1]
radius_outer = radius_inner;

% this is [ m, n]
% m means cell number, n means methods number
radius = [];

taLimStep = 120000/219; % 547.9425

% store the cell characteristic
cell_chr = [];

cell_enu = [];

real_radius = [];

for ii = 1:length(cell_data)
    
    [cell_chr(ii, 1), cell_chr(ii,2)] = cellChr(cell_data{ii, CELL_FREQ}, cell_data{ii,CELL_CLASS});
    
    cell_enu(ii, :) = cell_data{ii, LAT_LONG}; 
    
    real_radius(ii) = cell_data{ii, REAL_RADIUS};

    % radius process
    if custom_radius
        radius(ii, index_custom_radius) = str2double(cell_data{ii, CUSTOM_RADIUS});
    end
    
    if ta_radius
        [radius_inner(ii), radius(ii, index_ta_radius)] = taAlgorithm(str2double(cell_data{ii, CELL_TALIM}));
    end
    
    
    if oh_propagation_model
        radius(ii, index_oh_propagation_model) = ohPropagationAlgorithm( ...
            str2double(cell_data{ii, CELL_POWER}), str2double(cell_data{ii, CELL_ACCMIN}));
    end
end


% this is [m, 1]
if isd_radius
    radius(:, index_isd_radius) = isdAlgorithm(cell_enu, ISD_NumberOfClosestCells, Adjustment_Clearance);
end

if voronoi_radius_max || voronoi_radius_mean
    MAX_WAY = 1;
    MEAN_WAY = 2;
    
    v = VoronoiFortuneAlgo(cell_enu, 0.5);

    v.do();
    
    for ii = 1: length(cell_enu)
        p.x = cell_enu(ii, 1);
        p.y = cell_enu(ii, 2);
        
        if voronoi_radius_max
            radius(ii, index_voronoi_radius_max) = v.calculateRadius(p, MAX_WAY);
        end
        
        if voronoi_radius_mean
            radius(ii, index_voronoi_radius_mean) = v.calculateRadius(p, MEAN_WAY);
        end
    end
end

stat_report = zeros(1, size(radius, 2));
    
% after all the algorithm is processed.
% final calculate outer_radius & inner_radius
for ii = 1 : length(cell_data)
    
    %radius_outer(ii) = min([radius_custom(ii) radius_ta(ii) ...
    %    radius_oh(ii), radius_isd(ii)]);
    
    [radius_outer(ii), index] = min(radius(ii, :));
    
    stat_report(1, index) = stat_report(1, index) + 1;
    
    if radius_outer(ii)-2*taLimStep<radius_inner(ii)
        radius_inner(ii) = radius_outer(ii) - radius_outer(ii)-2.*taLimStep;
        if radius_inner(ii)<0
            radius_inner(ii) = 0;
        end
    end
end

% selected radius
methods = methods + 1;
index_selected_radius = methods;
radius(:, index_selected_radius) = radius_outer;
methods_name{index_selected_radius} = 'seleted radius';

% the real readius
methods = methods + 1;
index_tags = methods;
methods_name{index_tags} = 'real radius';

radius(:, index_tags) = real_radius;


% pie_legend = {'500', '1000', '3000', '10000', '20000', 'other'};
% M = [250, 750, 1250, 4750, 15250, 24750];
pie_legend = {'500', '3000', '10000', 'other'};
M = [250, 750, 5250, 14750];
for ii = 1 : size(radius, 2)-1
    relative_radius(:, ii) = radius(:, ii) - real_radius';
    
    figure;
    plot(relative_radius(:, ii));
    legend(methods_name{ii});
    
    figure;
    
    hist(abs(relative_radius(:, ii)));
    legend(methods_name{ii});
    
    figure;
    % 0~500 500~1000 1000~3000 3000~10000 10000~20000  20000~
    N = hist(abs(relative_radius(:, ii)), M);
 
    relative_count(:, ii) = N;
    
    pie(N);
    z = N==0;
    
    pie_index = 1;
    for jj = 1: length(z)
        if ~z(jj)
            legend_str{pie_index} = pie_legend{jj};
            pie_index = pie_index + 1;
        end
    end
    legend(legend_str, 'Location', 'SouthEast');
    title(methods_name{ii});

end

clear cell;
radius_by_cell = cell(5,3);
for ii = 1 : length(cell_chr)
    row_count = size(radius_by_cell{cell_chr(ii,1), cell_chr(ii,2)}, 1);
    radius_by_cell{cell_chr(ii,1), cell_chr(ii,2)}(row_count+1, :) = radius(ii, :);
end

relative_radius_by_cell = cell(5,3);
for ii = 1 : size(radius_by_cell, 1)
    for jj = 1 : size(radius_by_cell, 2)
        col_count = size(radius_by_cell{ii, jj}, 2);
        
        if col_count>0 
            for kk = 1 : col_count-1
                relative_radius_by_cell{ii, jj}(:, kk) = radius_by_cell{ii,jj}(:, kk) ...
                   - radius_by_cell{ii, jj}(:, col_count); 
            end
        end
    end
end

% caluculate the similar value
similar = [size(radius, 2), 1];
for ii = 1 : size(radius, 2)-1
    similar(ii) = similarAlgorithm(real_radius', radius(:, ii));
end

figure;
plot(similar, 'o--');
for ii = 1: length(similar)
    text(ii+0.1, similar(ii), methods_name{ii});
end
axis([0, length(similar)+1, 0, 1]);
grid on;
% legend show;


% output the stat info
figure;
bar(stat_report);

% output cell info
% figure;
% 
% plot(relative_radius(:, index_ta_radius:index_oh_propagation_model));
% legend(methods_name{index_ta_radius}, methods_name{index_oh_propagation_model});

figure;
plot(relative_radius(:, index_isd_radius:index_voronoi_radius_mean));
legend(methods_name{index_isd_radius}, methods_name{index_voronoi_radius_max}, methods_name{index_voronoi_radius_mean});


% figure;
% hold on;
% for ii = 1 : length(cell_enu)
% %ii = 2;
%     % draw start line
%     start_angle = degtorad(90-cell_angle(ii, 1));
%     stop_angle = degtorad(90-cell_angle(ii, 2));
%     
%     if start_angle<0
%         start_angle = start_angle + 2.*pi;
%     end
%     
%     if stop_angle<0
%         stop_angle = stop_angle + 2.*pi;
%     end
%     
%     if start_angle>=stop_angle
%         stop_angle = stop_angle + 2.*pi;
%     end
%     
%     x = linspace(cell_enu(ii, 1), cell_enu(ii, 1) + radius_outer(ii)*cos(start_angle), 1000);
%     y = tan(start_angle)*(x-cell_enu(ii, 1)) + cell_enu(ii, 2);
%     plot(x, y);
%     
%     % draw the stop line
%     x = linspace(cell_enu(ii, 1), cell_enu(ii, 1) + radius_outer(ii)*cos(stop_angle), 1000);
%     y = tan(stop_angle)*(x-cell_enu(ii, 1)) + cell_enu(ii, 2);
%     plot(x, y);
%     
%     % draw arc - outer radius
%     angle = linspace(start_angle, stop_angle, 1000);
%     plot(cell_enu(ii, 1) + radius_outer(ii).*cos(angle), cell_enu(ii, 2) + radius_outer(ii).*sin(angle));
%     
%     % draw arc - inner radius
%     angle = linspace(start_angle, stop_angle, 1000);
%     plot(cell_enu(ii, 1) + radius_inner(ii).*cos(angle), cell_enu(ii, 2) + radius_inner(ii).*sin(angle));
%     
% end
% hold off;
% axis equal;

