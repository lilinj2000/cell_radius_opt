function cell_data = initialCellData()
% inital cell data
% 

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


%  the fields of the cell_data as below:
LAC_CI = 1;
LAT_LONG = 2;
CELL_ANGLE = 3;
CELL_FREQ = 4;
CELL_CLASS = 5;
CELL_POWER = 6;
CUSTOM_RADIUS = 7;
CELL_TALIM = 8;
CELL_ACCMIN = 9;
REAL_RADIUS = 10;

MAX_FIELD = 10;

% gsm_file = 'GSM.cnai';
% gsm_file = 'all.cnai';

rogers = 1;
turkey = 2;
india = 3;

country = rogers;


switch country
    case {rogers}
        load rogers_gsm; % it's standard header
        % Header is below:
        % 1~5   BSC	SITE	altitude	latitude	longitude
        % 6~10  CELL	accmin	antenna_gain	antenna_tilt	antenna_type
        % 11~15 bcc	bcchno	bspwr	bspwrb	c_sys_type	
        % 16~20 cell_dir	cell_type	ci	env_char	height	
        % 21~25 lac latitude	longitude	max_altitude	max_cell_radius	
        % 26~30 mcc min_altitude	mnc	ncc	sector_angle	
        % 31~32 talim override_ta

        % Header enu
        LATITUDE = 22;
        LONGITUDE = 23;
        ACCMIN = 7;
        BCCHNO = 12;
        BSPWR = 13;
        C_SYS_TYPE = 15;
        CELL_DIR = 16;
        CELL_TYPE = 17;
        ANTENNA_TYPE = 10;
        MAX_CELL_RADIUS = 25;
        SECTOR_ANGLE = 30;
        TALIM = 31;

        MCC = 26;
        MNC = 28;
        LAC = 21;
        CI = 18;
    case {turkey}
        load turkey_gsm;
        % 1~5  BSC	SITE	latitude	longitude	CELL	
        % 6~10 accmin	antenna_type	bspwr	cell_dir	ci	
        % 11~15 lac	mcc	mnc	talim	sector_angle	
        % 16~20 max_cell_radius	ncc	bcc	bcch	latitude	
        % 21 longitude
        LATITUDE = 20;
        LONGITUDE = 21;
        ACCMIN = 6;
        BSPWR = 8;
        CELL_DIR = 9;
        % CELL_TYPE = 17;
        ANTENNA_TYPE = 7;
        MAX_CELL_RADIUS = 16;
        SECTOR_ANGLE = 15;
        TALIM = 14;

        MCC = 12;
        MNC = 13;
        LAC = 11;
        CI = 10;
    case {india}
        load india_gsm;
        
        % 1~5   CELL	LAT	LON	ANT_DIRECTION	ANT_BEAM_WIDTH
        % 6~10  ARFCN	BSIC	MCC	MNC	LAC
        % 11 CI
        LATITUDE = 2;
        LONGITUDE = 3;
        %ACCMIN = 6;
        %BSPWR = 8;
        %CELL_DIR = 9;
        %CELL_TYPE = 17;
        %ANTENNA_TYPE = 10;
        %MAX_CELL_RADIUS = 25;
        %SECTOR_ANGLE = 30;
        %TALIM = 31;
        
        MCC = 8;
        MNC = 9;
        LAC = 10;
        CI = 11;
end

% gsm_data = ImportGSMCellData(gsm_file);

%%%%%%%%%%%%%%Rogers Data%%%%%%%%%%%%%%%%%%%
% rural orange ville
% gsm_data = rural_orangeville_cell;
% gsm_tags = rural_orangeville_tags;

% suburban brampton 
% gsm_data = suburban_brampton_cell;
% gsm_tags = suburban_brampton_tags;

% suburban scarborough
% gsm_data = suburban_scarborough_cell;
% gsm_tags = suburban_scarborough_tags;

% suburan brampton + scarborough
% gsm_data = [suburban_brampton_cell'; suburban_scarborough_cell'];
% gsm_tags = [suburban_brampton_tags; suburban_scarborough_tags];

% urban dt
gsm_data = urban_dt_cell;
gsm_tags = urban_dt_tags;


% urban + suburban + rural
% gsm_data = [urban_dt_cell'; suburban_brampton_cell'; ...
%             suburban_scarborough_cell'; rural_orangeville_cell'];
% gsm_tags = [ urban_dt_tags; suburban_brampton_tags; ...
%             suburban_scarborough_tags; rural_orangeville_tags];


%%%%%%%%%%%%Turkey Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% turkey beyoglu
% gsm_data = turkey_beyoglu_cell;
% gsm_tags = turkey_beyoglu_tags2;

% turkey incek
% gsm_data = turkey_incek_cell;
% gsm_tags = turkey_incek_tags2;
% 
% % turkey sariyer
% gsm_data = turkey_sariyer_cell;
% gsm_tags = turkey_sariyer_tags2;
% 
% % turkey ulus
% gsm_data = turkey_ulus_cell;
% gsm_tags = turkey_ulus_tags;
% 
% % turkey uskudar
% gsm_data = turkey_uskudar_cell;
% gsm_tags = turkey_uskudar_tags2;
% 
% % turkey all
% gsm_data = [turkey_beyoglu_cell'; turkey_incek_cell'; ...
%             turkey_sariyer_cell'; turkey_ulus_cell'; ...
%             turkey_uskudar_cell'];
% gsm_tags = [ turkey_beyoglu_tags2; turkey_incek_tags2; ...
%             turkey_sariyer_tags2; turkey_ulus_tags; ...
%             turkey_uskudar_tags2];


%%%%%%%%%%%%India Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% india medium
% gsm_data = india_medium_cell;
% gsm_tags = india_medium_tags;

% india rural
% gsm_data = india_rural_cell;
% gsm_tags = india_rural_tags2;

% india rural highway
% gsm_data = india_rural_highway_cell;
% gsm_tags = india_rural_highway_tags2;

% india all
% gsm_data = [india_medium_cell'; india_rural_cell'; ...
%             india_rural_highway_cell'];
% gsm_tags = [ india_medium_tags; india_rural_tags2; ...
%             india_rural_highway_tags2];



% cell type 'MACRO', 'MICRO', 'PICO'
% cell_type = 'MACRO';
% index_cell = 0;
% for ii = 1 : length(gsm_data)
%     if strcmp(gsm_data{ii}(CELL_TYPE), cell_type)
%         index_cell = index_cell + 1;
%         gsm_data_internal{index_cell} = gsm_data{ii};
%     end
% end
% gsm_data = gsm_data_internal;

clear cell;
cell_data = cell(0, MAX_FIELD);
tags = [length(gsm_tags), 1];

for ii = 1 : length(gsm_data)
    switch country
        case {rogers, turkey}
            
            [cell_angle(1, 1), cell_angle(1, 2)] = sectorAngle(gsm_data{ii}{ANTENNA_TYPE}, ...
            str2double(gsm_data{ii}{CELL_DIR}), str2double(gsm_data{ii}{SECTOR_ANGLE}));
        
            cell_freq = gsm_data{ii}(C_SYS_TYPE);
            cell_type = gsm_data{ii}(CELL_TYPE);
            cell_power = gsm_data{ii}{BSPWR};
            custom_radius = gsm_data{ii}{MAX_CELL_RADIUS};
            cell_talim = gsm_data{ii}{TALIM};
            cell_accmin = gsm_data{ii}{ACCMIN};
            
            
        case {india}
            
            cell_angle = [];
            
            
            
    end
    
    lac_ci = strcat(gsm_data{ii}{LAC}, '-', ...
                gsm_data{ii}{CI});
            
    cell_latlong(1, 1) = dms2deg(gsm_data{ii}{LATITUDE});
    cell_latlong(1, 2) = dms2deg(gsm_data{ii}{LONGITUDE});

    cell_lat_long_rad = degtorad(cell_latlong);
    if ii==1
        cell_lat_long_rad_base = cell_lat_long_rad;
    end
    cell_enu = convertlatlong2enu(cell_lat_long_rad, cell_lat_long_rad_base);
    
    cell_data{ii, LAC_CI} = lac_ci;
    cell_data{ii, LAT_LONG} = cell_enu;
    cell_data{ii, CELL_ANGLE} = cell_angle;
    cell_data{ii, CELL_FREQ} = cell_freq;
    cell_data{ii, CELL_CLASS} = cell_type;
    cell_data{ii, CELL_POWER} = cell_power;
    cell_data{ii, CUSTOM_RADIUS} = custom_radius;
    cell_data{ii, CELL_TALIM} = cell_talim;
    cell_data{ii, CELL_ACCMIN} = cell_accmin;
    
    
    cell1 = strcat(gsm_data{ii}{MCC}, '-', ...
                gsm_data{ii}{MNC}, '-', ...
                gsm_data{ii}{LAC}, '-', ...
                gsm_data{ii}{CI});
    
    for jj = 1 : length(gsm_tags);
        cell2 = gsm_tags{jj, 2};
        if strcmp(cell1, cell2)
            tags(ii) = gsm_tags{jj, 5}(4);
            break;
        end
    end
    
    cell_data{ii, REAL_RADIUS} = tags(ii);
    
end

%check tags
if ~isempty(find(tags==0))
    error('Wrong tags!!!');
end

end