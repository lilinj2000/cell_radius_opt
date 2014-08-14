clc;
% clear;

% fWriteId = fopen('Rural_Orangeville.txt', 'w+');

% india_medium_tags = cell();
% index = 1;
% for ii = 1: length(Tags_uskudar)
% 
%     if isempty(strfind(Tags_uskudar{ii,2}, ','))
%         turkey_uskudar_tags(index, :) = Tags_uskudar(ii, :);
%         index = index + 1;
%     end
% end

% fclose(fWriteId);

% gsm_file = 'celldata_Rogers.cnai';
% rogers_cell_data = ImportGSMCellData('celldata_Rogers.cnai');

% fReadId = fopen('Rural_Orangeville.txt', 'r');
% 
% tag_data = textscan(fReadId, '%s %d');
% 
% fclose(fReadId);

% load rogers_cell_data;
% load urban_dt_tags;
% load suburban_brampton_tags;
% load suburban_scarborough_tags;
% load rural_orangeville_tags;

% Header is below:
% 1~5   BSC	SITE	altitude	latitude	longitude
% 6~10  CELL	accmin	antenna_gain	antenna_tilt	antenna_type
% 11~15 bcc	bcchno	bspwr	bspwrb	c_sys_type	
% 16~20 cell_dir	cell_type	ci	env_char	height	
% 21~25 lac latitude	longitude	max_altitude	max_cell_radius	
% 26~30 mcc min_altitude	mnc	ncc	sector_angle	
% 31~32 talim override_ta

% MCC = 26;
% MNC = 28;
% LAC = 21;
% CI = 18;
% 
%    
% index_urban = 1;
% index_suburban_scarborough = 1;
% index_suburban_brampton = 1;
% index_rural_orangeville = 1;



% for jj = 1 : length(rogers_cell_data)
%     cell1 = strcat(rogers_cell_data{jj}(MCC), '-', ...
%         rogers_cell_data{jj}(MNC), '-', ...
%         rogers_cell_data{jj}(LAC), '-', ...
%         rogers_cell_data{jj}(CI));
%     
%     cells(jj) = cell1;
% end

% 1~5   CELL	LAT	LON	ANT_DIRECTION	ANT_BEAM_WIDTH
% 6~10  ARFCN	BSIC	MCC	MNC	LAC
% 11 CI
% MCC = 8;
% MNC = 9;
% LAC = 10;
% CI = 11;


% 1~5  BSC	SITE	latitude	longitude	CELL	
% 6~10 accmin	antenna_type	bspwr	cell_dir	ci	
% 11~15 lac	mcc	mnc	talim	sector_angle	
% 16~20 max_cell_radius	ncc	bcc	bcch	latitude	
% 21 longitude
MCC = 12;
MNC = 13;
LAC = 11;
CI = 10;

% index = 1;
% for jj = 1 : length(turkey_cell)
%     cell1 = strcat(turkey_cell{jj}{MCC}, '-', ...
%         turkey_cell{jj}{MNC}, '-', ...
%         turkey_cell{jj}{LAC}, '-', ...
%         turkey_cell{jj}{CI});
%     
% %     urban_dt_cell_found = 0;
%     for ii = 1: length(turkey_uskudar_tags)
%         cell2 = turkey_uskudar_tags{ii,2};
%         
%         if strcmp(cell1, cell2)
%             turkey_uskudar_cell(index) = turkey_cell(jj);
%             index = index + 1;
% %             urban_dt_cell_found = 1;
%             break;
%         end
%     end
% end

index = 1;
turkey_uskudar_tags2 = cell(0);
 for ii = 1: length(turkey_uskudar_tags)
    cell2 = turkey_uskudar_tags{ii,2};

    for jj = 1 : length(turkey_uskudar_cell)
        cell1 = strcat(turkey_uskudar_cell{jj}{MCC}, '-', ...
        turkey_uskudar_cell{jj}{MNC}, '-', ...
        turkey_uskudar_cell{jj}{LAC}, '-', ...
        turkey_uskudar_cell{jj}{CI});
    
        if strcmp(cell1, cell2)
            turkey_uskudar_tags2(index, :) = turkey_uskudar_tags(ii, :);
            index = index + 1;
    %             urban_dt_cell_found = 1;
            break;
        end

    end
    
end




% for ii = 1: length(turkey_uskudar_tags)
%     cell = turkey_uskudar_tags{ii,2};
%     
%     fields = regexp(cell, '-', 'split');
%     
%     cell = strcat('286-1-', fields{3}, '-', fields{4});
%     
%     turkey_uskudar_tags{ii, 2} = cell;
%     turkey_uskudar_tags{ii, 3} = cell;
%     
% end

% for jj = 1 : length(rogers_cell_data)
%     cell1 = strcat(rogers_cell_data{jj}(MCC), '-', ...
%         rogers_cell_data{jj}(MNC), '-', ...
%         rogers_cell_data{jj}(LAC), '-', ...
%         rogers_cell_data{jj}(CI));
%     
%     urban_dt_cell_found = 0;
%     for ii = 1: length(urban_dt_tags)
%         cell2 = urban_dt_tags{ii,2};
%         
%         if strcmp(cell1, cell2)
%             urban_dt_cell(index_urban) = rogers_cell_data(jj);
%             index_urban = index_urban + 1;
%             urban_dt_cell_found = 1;
%             break;
%         end
%     end
%     
%     if ~urban_dt_cell_found
%         suburban_scarborough_cell_found = 0;
%         for ii = 1: length(suburban_scarborough_tags)
%             cell2 = suburban_scarborough_tags{ii,2};
%         
%             if strcmp(cell1, cell2)
%                 suburban_scarborough_cell(index_suburban_scarborough) = rogers_cell_data(jj);
%                 index_suburban_scarborough = index_suburban_scarborough + 1;
%                 suburban_scarborough_cell_found = 1;
%                 break;
%             end
%         end
%     end
%     
%     if ~urban_dt_cell_found && ~suburban_scarborough_cell_found
%         suburban_brampton_cell_found = 0;
%         for ii = 1: length(suburban_brampton_tags)
%             cell2 = suburban_brampton_tags{ii,2};
%         
%             if strcmp(cell1, cell2)
%                 suburban_brampton_cell(index_suburban_brampton) = rogers_cell_data(jj);
%                 index_suburban_brampton = index_suburban_brampton + 1;
%                 suburban_brampton_cell_found = 1;
%                 break;
%             end
%         end
%     end
%     
%     if ~urban_dt_cell_found && ~suburban_scarborough_cell_found && ~suburban_brampton_cell_found
%         rural_orangeville_cell_found = 0;
%         for ii = 1: length(rural_orangeville_tags)
%             cell2 = rural_orangeville_tags{ii,2};
%         
%             if strcmp(cell1, cell2)
%                 rural_orangeville_cell(index_rural_orangeville) = rogers_cell_data(jj);
%                 index_rural_orangeville = index_rural_orangeville + 1;
%                 rural_orangeville_cell_found = 1;
%                 break;
%             end
%         end
%     end
% end


