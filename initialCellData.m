function cell_data = initialCellData(country, data_area)
% inital cell data
% 

importSysVar;

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
        HEIGHT = 20;
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
        
%     case {china}

end

% gsm_data = ImportGSMCellData(gsm_file);

switch data_area
%     case rogers_area_rural_orangeville
% 
%         gsm_data = rural_orangeville_cell;
%         gsm_tags = rural_orangeville_tags;
% 
%     case rogers_area_suburban_brampton
%  
%         gsm_data = suburban_brampton_cell;
%         gsm_tags = suburban_brampton_tags;
% 
%     case rogers_area_suburban_scarborough
% 
%         gsm_data = suburban_scarborough_cell;
%         gsm_tags = suburban_scarborough_tags;
% 
%     case rogers_area_suburban
%  
%         gsm_data = [suburban_brampton_cell'; suburban_scarborough_cell'];
%         gsm_tags = [suburban_brampton_tags; suburban_scarborough_tags];
% 
%     case rogers_area_urban_dt
% 
%         gsm_data = urban_dt_cell;
%         gsm_tags = urban_dt_tags;

    case rogers_area_all
        % urban + suburban + rural
%         gsm_data = [urban_dt_cell'; suburban_brampton_cell'; ...
%                     suburban_scarborough_cell'; rural_orangeville_cell'];
        gsm_data = rogers_cell_data;
        gsm_tags = [ urban_dt_tags; suburban_brampton_tags; ...
                    suburban_scarborough_tags; rural_orangeville_tags];
                
%     case turkey_beyoglu
% 
%     gsm_data = turkey_beyoglu_cell;
%     gsm_tags = turkey_beyoglu_tags2;
% 
%     case turkey_incek
%     gsm_data = turkey_incek_cell;
%     gsm_tags = turkey_incek_tags2;
% 
% 
%     case turkey_sariyer
% 
%     gsm_data = turkey_sariyer_cell;
%     gsm_tags = turkey_sariyer_tags2;
% 
%     case turkey_ulus
% 
%     gsm_data = turkey_ulus_cell;
%     gsm_tags = turkey_ulus_tags;
% 
%     case turkey_uskudar
% 
%     gsm_data = turkey_uskudar_cell;
%     gsm_tags = turkey_uskudar_tags2;
    
    case turkey_all
    % turkey all
%     gsm_data = [turkey_beyoglu_cell'; turkey_incek_cell'; ...
%                 turkey_sariyer_cell'; turkey_ulus_cell'; ...
%                 turkey_uskudar_cell'];
            
    gsm_data = turkey_cell;
    gsm_tags = [ turkey_beyoglu_tags2; turkey_incek_tags2; ...
                turkey_sariyer_tags2; turkey_ulus_tags; ...
                turkey_uskudar_tags2];

%     case india_medium
    %%%%%%%%%%%%India Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % india medium
%     gsm_data = india_medium_cell;
%     gsm_tags = india_medium_tags;
% 
%     case india_rural
%     % india rural
%     gsm_data = india_rural_cell;
%     gsm_tags = india_rural_tags2;
% 
%     case india_rural_highway
%     % india rural highway
%     gsm_data = india_rural_highway_cell;
%     gsm_tags = india_rural_highway_tags2;

    case india_all
    % india all
%     gsm_data = [india_medium_cell'; india_rural_cell'; ...
%                 india_rural_highway_cell'];
    gsm_data = india_cell';
    gsm_tags = [ india_medium_tags; india_rural_tags2; ...
                india_rural_highway_tags2];
            
    case china_beijing
        load beijing_wcdma;
        
        LATITUDE = 1;
        LONGITUDE = 2;
        CELL_ID = 14;
        
%         gsm_data = beijing_wcdma_cells;
        gsm_data = unq_beijing_cell;
        gsm_tags = beijing_wcdma_tags;
            
    case china_wuxi
        load wx_wcdma_cell;
        
        LATITUDE = 2;
        LONGITUDE = 3;
        CELL_ID = 1;
        
%         gsm_data = wx_wcdma_cells;
        gsm_data = wx_new_total_data;
        gsm_tags = wx_wcdma_tags;
        
    otherwise
        error('invalid input');
% suburban + rural
% gsm_data = [suburban_brampton_cell'; ...
%             suburban_scarborough_cell'; rural_orangeville_cell'];
% gsm_tags = [ suburban_brampton_tags; ...
%             suburban_scarborough_tags; rural_orangeville_tags];
end



cell_data = cell(0, LAST_FIELD);
tags = zeros(length(gsm_data), 1);
index_cell_data = 0;
for ii = 1 : length(gsm_data)
    switch country
        case {rogers, turkey}
            
            [cell_angle(1, 1), cell_angle(1, 2)] = sectorAngle(gsm_data{ii}{ANTENNA_TYPE}, ...
            str2double(gsm_data{ii}{CELL_DIR}), str2double(gsm_data{ii}{SECTOR_ANGLE}));
        
            if country==turkey
                cell_freq = [];
                cell_type = [];
                cell_height = [];
            else
                cell_freq = gsm_data{ii}{C_SYS_TYPE};
                cell_type = gsm_data{ii}{CELL_TYPE};
                cell_height = gsm_data{ii}{HEIGHT};
            end
            
            
            cell_power = gsm_data{ii}{BSPWR};
            custom_radius = gsm_data{ii}{MAX_CELL_RADIUS};
            cell_talim = gsm_data{ii}{TALIM};
            cell_accmin = gsm_data{ii}{ACCMIN};
            
            
            
        case {india, china}
            
            cell_angle = [];
            cell_freq = [];
            cell_type = [];
            cell_height = [];
            cell_power = [];
            custom_radius = [];
            cell_talim = [];
            cell_accmin = [];
            
            
            
    end
    
    if data_area==china_beijing
        lac = fix(gsm_data(ii,CELL_ID)/2^16);
        cid = mod(gsm_data(ii, CELL_ID), 2^16);
        
        lac_ci = strcat(num2str(lac), '-', num2str(cid));
        cell_latlong(1, 1) = gsm_data(ii,LATITUDE);
        cell_latlong(1, 2) = gsm_data(ii,LONGITUDE);
    elseif data_area==china_wuxi
        
        lac_ci = strcat('0-', num2str(gsm_data(ii, CELL_ID)));
        cell_latlong(1, 1) = gsm_data(ii,LATITUDE);
        cell_latlong(1, 2) = gsm_data(ii,LONGITUDE);
    else
        lac_ci = strcat(gsm_data{ii}{LAC}, '-', ...
                gsm_data{ii}{CI});
            
        latitude = gsm_data{ii}{LATITUDE};
        longitude = gsm_data{ii}{LONGITUDE};
        
        if country~=india
            if isempty(latitude) || isempty(longitude) ...
                    || length(latitude)<=4 || length(longitude)<=4
                continue;
            end
        end
        
        cell_latlong(1, 1) = dms2deg(gsm_data{ii}{LATITUDE});
        cell_latlong(1, 2) = dms2deg(gsm_data{ii}{LONGITUDE});
    end
            
    

    cell_lat_long_rad = degtorad(cell_latlong);
    if ii==1
        cell_lat_long_rad_base = cell_lat_long_rad;
    end
    cell_enu = convertlatlong2enu(cell_lat_long_rad, cell_lat_long_rad_base);
    
    
    index_cell_data =  index_cell_data + 1;
    
    cell_data{index_cell_data, LAC_CI} = lac_ci;
    cell_data{index_cell_data, LAT_LONG} = cell_enu;
    cell_data{index_cell_data, CELL_ANGLE} = cell_angle;
    cell_data{index_cell_data, CELL_FREQ} = cell_freq;
    cell_data{index_cell_data, CELL_CLASS} = cell_type;
    cell_data{index_cell_data, CELL_POWER} = cell_power;
    cell_data{index_cell_data, CUSTOM_RADIUS} = custom_radius;
    cell_data{index_cell_data, CELL_TALIM} = cell_talim;
    cell_data{index_cell_data, CELL_HEIGHT} = cell_height;
    cell_data{index_cell_data, CELL_ACCMIN} = cell_accmin;
    
    
    cell1 = lac_ci;
    
    for jj = 1 : length(gsm_tags);
        cell2_string = gsm_tags{jj, 2};
        cell2_split = regexp(cell2_string, '-', 'split');
        
        cell2 = strcat(cell2_split{3}, '-', cell2_split{4});
        
        if strcmp(cell1, cell2)
            
            tags(ii) = gsm_tags{jj, 5}(4);
            
            break;
        end
    end
    
    cell_data{index_cell_data, REAL_RADIUS} = tags(ii);
    
end

% remove the invalid tags
% invalid_tags = find(tags==0);
% cell_data(invalid_tags, :) = [];

% if invalid_(find(tags==0))
%     error('Wrong tags!!!');
% end

end