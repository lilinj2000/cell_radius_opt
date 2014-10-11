%% CELLRADIUSPROCESS
% the cell radius optimization main process


%% environment initial
% close all;
% clear;
% clc;


%% Parameter Confgure


%%% General parametr configure

% Specifies the distance between one site and another when searching 
% for co-located cells. 
% If the distance is less than the value of Adjustment_Clearance, 
% the sites are to be converged into one site.
% Unit: meters
Adjustment_Clearance = 50; 

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

%%% ISD & Voronoi Algorithm configure

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

%% Algorithm Method Selection

custom_method = false;
ta_method = false;
oh_method = false;
toh_method = true;
isd_method = false;
voronoi_method = false;

methods = 0;
index_custom_radius = 0;
index_ta_radius = 0;
index_oh_model = 0;
index_toh_model = 0;
index_isd_max_radius = 0;
index_isd_mean_radius = 0;
index_isd_min_radius = 0;
index_isd_angle_radius = 0;

index_voronoi_vradius_max = 0;
index_voronoi_vradius_mean = 0;
index_voronoi_vradius_min = 0;
index_voronoi_vradius_angle = 0;
index_voronoi_sradius_max = 0;
index_voronoi_sradius_mean = 0;
index_voronoi_sradius_min = 0;
index_voronoi_sradius_angle = 0;


methods_name = {};

if custom_method
    methods = methods + 1;
    index_custom_radius = methods;
    methods_name{index_custom_radius} = 'custom radius';
end

if ta_method 
    methods = methods + 1;
    index_ta_radius = methods;
    methods_name{index_ta_radius} = 'TA';
end

if oh_method
    methods = methods + 1;
    index_oh_model = methods;
    methods_name{index_oh_model} = 'OH'; 
end

if toh_method
    methods = methods + 1;
    index_toh_model = methods;
    methods_name{index_toh_model} = 'TOH';
end

if isd_method
    methods = methods + 1;
    index_isd_max_radius = methods;
    methods_name{index_isd_max_radius} = 'ISD max';
    
    methods = methods + 1;
    index_isd_mean_radius = methods;
    methods_name{index_isd_mean_radius} = 'ISD mean';
    
    methods = methods + 1;
    index_isd_min_radius = methods;
    methods_name{index_isd_min_radius} = 'ISD min';
    
    methods = methods + 1;
    index_isd_angle_radius = methods;
    methods_name{index_isd_angle_radius} = 'ISD angle';
end

if voronoi_method
    methods = methods + 1;
    index_voronoi_vradius_max = methods;
    methods_name{index_voronoi_vradius_max} = 'VORONOI vmax';
    
    methods = methods + 1;
    index_voronoi_vradius_mean = methods;
    methods_name{index_voronoi_vradius_mean} = 'VORONOI vmean';
    
    methods = methods + 1;
    index_voronoi_vradius_min = methods;
    methods_name{index_voronoi_vradius_min} = 'VORONOI vmin';
    
    methods = methods + 1;
    index_voronoi_vradius_angle = methods;
    methods_name{index_voronoi_vradius_angle} = 'VORONOI vangle';
    
    methods = methods + 1;
    index_voronoi_sradius_max = methods;
    methods_name{index_voronoi_sradius_max} = 'VORONOI smax';
    
    methods = methods + 1;
    index_voronoi_sradius_mean = methods;
    methods_name{index_voronoi_sradius_mean} = 'VORONOI smean';
    
    methods = methods + 1;
    index_voronoi_sradius_min = methods;
    methods_name{index_voronoi_sradius_min} = 'VORONOI smin';
    
    methods = methods + 1;
    index_voronoi_sradius_angle = methods;
    methods_name{index_voronoi_sradius_angle} = 'VORONOI sangle';
end

%% include the custom variables
importSysVar;

% urban = 1;
% suburban = 2;
% rural = 3;
% all = 4;
% area_type = all;

subhead = {'(turkey all)'};

skip_power = true;
skip_angle = true;

% initial the cell data
% cell_data = initialCellData(india, india_all);
% load rogers_radius;
% load turkey_radius;
% load beijing_radius;

radius_ready = false;

idx_real_radius = 4;

if radius_ready
    
%     radius = [rogers_ta_radius, rogers_oh_radius, rogers_toh_radius, rogers_isd_radius, rogers_voro_radius];
    
%     radius = [urban_dt_ta_radius, urban_dt_oh_radius, urban_dt_toh_radius, urban_dt_isd_radius, urban_dt_voro_radius];
%     real_radius = urban_dt_radius(:, idx_real_radius);

    %     radius = [wx_isd_radius, wx_voro_radius];

%     radius = [suburban_brampton_ta_radius, suburban_brampton_oh_radius, suburban_brampton_toh_radius, suburban_brampton_isd_radius, suburban_brampton_voro_radius];
%     radius = [radius; suburban_scarborough_ta_radius, suburban_scarborough_oh_radius, suburban_scarborough_toh_radius, suburban_scarborough_isd_radius, suburban_scarborough_voro_radius];
%     real_radius = suburban_brampton_radius(:, idx_real_radius);
%     real_radius = [real_radius; suburban_scarborough_radius(:,idx_real_radius)];

%     radius = [rural_orangeville_ta_radius, rural_orangeville_oh_radius, rural_orangeville_toh_radius, rural_orangeville_isd_radius, rural_orangeville_voro_radius];
%     real_radius = rural_orangeville_radius(:, idx_real_radius);
    
% radius = [turkey_ta_radius, turkey_oh_radius, turkey_isd_radius, turkey_voro_radius];

radius = [beijing_isd_radius, beijing_voro_radius];
radius = radius.*3.6 + 7500;
    %     radius = rogers_isd_radius;

%     radius = [wx_isd_radius, wx_voro_radius];

    real_radius = r_real_radius(:, idx_real_radius);
    drawAnalysisResult(radius, real_radius, methods_name, subhead);
    return;
end

cell_data = turkey_cell;
cell_data_real_radius = turkey_cell_info;



% [ m, 1]
radius_inner = zeros(length(cell_data), 1);

% the real radius is [ m, 1]
radius_outer = radius_inner;

% this is [ m, n]
% m means cell number, n means methods number
radius = [];

% taLimStep = 120000/219; % 547.9425

% store the cell characteristic
% t_lac_ci = cell(0);
% t_cell_enu = [];
% 
% t_real_radius = [];
% t_cell_power = [];
% t_cell_angle = [];
% 
% r_lac_ci = cell(0);
% r_cell_enu = [];
% r_real_radius = [];
% r_cell_power = [];
% r_cell_angle = [];
% 
% 
% for ii = 1:length(cell_data)
%     
%     t_lac_ci(ii, 1) = cellstr(cell_data{ii, LAC_CI});
%     t_cell_enu(ii, :) = cell_data{ii, LAT_LONG}; 
%     t_real_radius(ii, 1) = cell_data{ii, REAL_RADIUS};
%     
%     if ~isempty(cell_data_real_radius{ii, 3})
%         t_real_radius(ii, [2:4]) = cell_data_real_radius{ii, 3};
%     end
%     
%     if ~skip_power
%         t_cell_power(ii, 1) = str2double(cell_data{ii, CELL_POWER});
%         t_cell_power(ii, 2) = str2double(cell_data{ii, CELL_ACCMIN});
%     end
%    
%     
%    if ~skip_angle
%         t_cell_angle(ii, :) = cell_data{ii, CELL_ANGLE};
%    end
%     
%     
%     if t_real_radius(ii, 1)~=0
%         
%         r_lac_ci = [r_lac_ci; cellstr(t_lac_ci(ii, :))];
%         r_cell_enu = [r_cell_enu; t_cell_enu(ii, :)];
%         r_real_radius = [r_real_radius; t_real_radius(ii, :)];
%         
%         if ~skip_power
%             r_cell_power = [r_cell_power; t_cell_power(ii, :)];
%         end
%         
%         if ~skip_angle
%             r_cell_angle = [r_cell_angle; t_cell_angle(ii, :)];
%         end
%         
%         idx_radius = length(r_real_radius);
%         if custom_method
%             cust_radius = str2double(cell_data{ii, CUSTOM_RADIUS});
%             
%             radius(idx_radius, index_custom_radius) =  cust_radius;
%         end
%         
%         if ta_method
%             
%             [ta_radius_inner, ta_radius_outer] = taAlgorithm(str2double(cell_data{ii, CELL_TALIM}));
%             
%             radius_inner(idx_radius) = ta_radius_inner;
%             radius(idx_radius, index_ta_radius) = ta_radius_outer;
%         end
%     end
%     
%     
%     %[enu_lat, enu_long, start_angle, stop_angle, bspwr, accmin]
%         
% end



% toh_model = [0.1043, 150.1801];
% toh_model =
% 
%     0.1043  150.1801

% toh_model = [1.8964, 146.9522];
% toh_model =
% 
%     1.8964  146.9522
    
% toh_model = [];
if oh_method
    
    oh_cell_power = r_cell_power;
    radius(:, index_oh_model) = ohPropagationAlgorithm(oh_cell_power(:, 1), oh_cell_power(:, 2));
    
    
end

if toh_method
    toh_cell_radius_power = [r_cell_power, r_real_radius(:, idx_real_radius)];
    [A, B, radius(:, index_toh_model)] = tohPropagationAlgorithm(toh_cell_radius_power);
    
    toh_model = [A, B];
end

% this is [m, 1]
if isd_method
%     isd_cell_enu_angle = [t_cell_enu, t_cell_angle, t_cell_power, t_real_radius];
    
    radiusISD = isdAlgorithm(t_cell_enu, t_cell_angle, t_cell_power, t_real_radius(:, idx_real_radius), ...
        ISD_NumberOfClosestCells, Adjustment_Clearance, toh_model);
    
    radius(:, index_isd_max_radius:index_isd_angle_radius) = radiusISD;
end

% if isd_max_radius
%     radius(:, index_isd_max_radius) = isdAlgorithm(cell_enu, ISD_NumberOfClosestCells, Adjustment_Clearance, ISD_MAX_WAY);
% end

if voronoi_method
        
%     load india_voronoi;
%     
%     v = india_voro;
% load rogers_voronoi;
% v = rogers_voro;
load turkey_voronoi;
v = turkey_voro;

%     v = VoronoiFortuneAlgo(t_cell_enu, 0.5)
%     v.do();
%     
%     v.splitSegList();
    
%     v_cell_enu_angle = [t_cell_enu, t_cell_angle, t_cell_power, t_real_radius];
    radius(:, index_voronoi_vradius_max:index_voronoi_sradius_angle) = v.calculateRadius(t_cell_enu, t_cell_angle, t_cell_power, t_real_radius(:, idx_real_radius), Adjustment_Clearance, toh_model);
    
end



% stat_report = zeros(1, size(radius, 2));
%     
% % after all the algorithm is processed.
% % final calculate outer_radius & inner_radius
% for ii = 1 : length(cell_data)
%     
%     %radius_outer(ii) = min([radius_custom(ii) radius_ta(ii) ...
%     %    radius_oh(ii), radius_isd(ii)]);
%     
%     [radius_outer(ii), index] = min(radius(ii, :));
%     
%     stat_report(1, index) = stat_report(1, index) + 1;
%     
%     if radius_outer(ii)-2*taLimStep<radius_inner(ii)
%         radius_inner(ii) = radius_outer(ii) - radius_outer(ii)-2.*taLimStep;
%         if radius_inner(ii)<0
%             radius_inner(ii) = 0;
%         end
%     end
% end

% selected radius
% methods = methods + 1;
% index_selected_radius = methods;
% radius(:, index_selected_radius) = radius_outer;
% methods_name{index_selected_radius} = 'seleted radius';

% remove the real radius is 0 data.
% [rows, cols] = find(real_radius'==0);
% real_radius(rows) = [];
% radius(rows, :) = [];
% 
% % remove the radius, according the area_type
% min_radius = 50;
% max_radius_urban = 1000;
% max_radius_suburban = 5000;
% max_radius_rural = 30000;
% 
% tmp_real_radius = [];
% tmp_radius = [];
% tmp_cell_power = [];
% for ii=1:length(real_radius)
% 
%     match = false;
%     switch area_type
%         case urban
%             if real_radius(ii)<=max_radius_urban ...
%                 && real_radius(ii)>=min_radius ...
%                 && radius(ii, eliminat_by_method) <= max_radius_urban
% 
%                 match = true;
%             end
% 
%         case suburban
%             if real_radius(ii)<=max_radius_suburban  ...
%                 && real_radius(ii)>max_radius_urban ...
%                 && radius(ii, eliminat_by_method) <= max_radius_suburban ...
%                 && radius(ii, eliminat_by_method) > max_radius_urban
% 
%                 match = true;
%             end
%         case rural
%             if real_radius(ii)<=max_radius_rural  ...
%                 && real_radius(ii)>max_radius_suburban ...
%                 && radius(ii, eliminat_by_method) > max_radius_suburban
% 
%                 match = true;
%             end
%             
%         case all
%             
%             if radius(ii, eliminat_by_method)<=max_radius_urban
%                
%                 if real_radius(ii)<=max_radius_urban ...
%                         && real_radius(ii)>=min_radius
%                     match = true;
%                 end
%                 
%             elseif radius(ii, eliminat_by_method)<=max_radius_suburban
%                 
%                 if real_radius(ii)>max_radius_suburban ...
%                     || real_radius(ii)<max_radius_urban
%                     match = false;
%                 else
%                     match = true;
%                 end
%                 
%             elseif radius(ii, eliminat_by_method)>max_radius_suburban ...
%                 
%                 if real_radius(ii)<max_radius_suburban ...
%                     || real_radius(ii)>max_radius_rural
%                     match = false;
%                 else
%                     match = true;
%                 end
%             end
%             
%             
%     end
% 
%     if match || (~eliminate_invalid_data)
%         index = length(tmp_real_radius) + 1;
%         tmp_real_radius(index) = real_radius(ii);
%         tmp_radius(index, :) = radius(ii, :);
%         
%         tmp_cell_power(index, :) = cell_power(ii, :);
%     end
% end
%     
% real_radius = tmp_real_radius;
% radius = tmp_radius;
% cell_power = tmp_cell_power;

% X = radius(:, index_oh_model:index_voronoi_sradius_angle);
% % power(10, (str2double(cell_data{ii, CELL_POWER}) ...
% %         + str2double(cell_data{ii, CELL_ACCMIN})+21-116)/35)
% % X(:, size(X, 2)+1) = power(10, (cell_power(:,1)+cell_power(:, 2)+21-116)/35);
% 
% [b,se,pval,inmodel,stats,nextstep,history] = stepwisefit(X, r_real_radius, 'display', 'off');
% 
% if find(inmodel>0)
%     methods = methods + 1;
%     index_regression = methods;
%     radius(:, index_regression) = stats.intercept;
%     methods_name{index_regression} = 'Regression';
%     for ii=1:length(inmodel)
% 
%         if inmodel(ii)~=0
%             radius(:, index_regression) = radius(:, index_regression) + b(ii).*X(:, ii);
%         end
%     end
% end

% if area_type==all %&& ~eliminate_invalid_data
%     methods = methods + 1;
%     index_best_method = methods;
%     methods_name{index_best_method} = 'Proposal';
%     % for rogers, the urban&suburban just use ISD radius,
%     % the rural, use Voronoi mean radius
%     for ii=1:length(radius)
%         if radius(ii, eliminat_by_method)<=max_radius_urban
%             radius(ii, index_best_method) = radius(ii, index_isd_mean_radius);
% %         elseif radius(ii, eliminat_by_method)<=max_radius_suburban
% %         if radius(ii, eliminat_by_method)<=max_radius_suburban
% %             radius(ii, index_best_method) = radius(ii, index_voronoi_radius_isd_max);
%         else
%             radius(ii, index_best_method) = radius(ii, index_isd_max_radius);
% %             radius(ii, index_best_method) = radius(ii, index_voronoi_radius_mean); 
% %             radius(ii, index_best_method) = 10320 + 0.87*radius(ii, index_voronoi_radius_min); 
%         end
%     end
% end





