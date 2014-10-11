% celld data process

%rogers
% LAC = 21;
% CI = 18;
% LATITUDE = 22;
% LONGITUDE = 23;
% gsm_data = rogers_cell_data;


%turkey
% LAC = 11;
% CI = 10;
% LATITUDE = 20;
% LONGITUDE = 21;
% gsm_data = turkey_data;

% india
% LATITUDE = 2;
% LONGITUDE = 3;
% LAC = 10;
% CI = 11;
% gsm_data = india_data;

% beijing
% LATITUDE = 1;
% LONGITUDE = 2;
% CELL_ID = 14;
% gsm_data = unq_beijing_cell;

%% wuxi
% LATITUDE = 3;
% LONGITUDE = 2;
% CELL_ID = 1;
% gsm_data = wx_new_total_data;
% 
% index = 0;
% cell_data = cell(0);
% for ii = 1:length(gsm_data)
%     
%     % beijing start
% %     lac = fix(gsm_data(ii,CELL_ID)/2^16);
% %     cid = mod(gsm_data(ii, CELL_ID), 2^16);
% % 
% %     lac_ci = strcat(num2str(lac), '-', num2str(cid));
% %     
% %     latitude = gsm_data(ii, LATITUDE);
% %     longitude = gsm_data(ii, LONGITUDE);
%     % beijing end
%     
%     % wuxi start
%     lac_ci = strcat('0-', num2str(gsm_data(ii, CELL_ID)));
%     latitude = gsm_data(ii, LATITUDE);
%     longitude = gsm_data(ii, LONGITUDE);
%     % wuxi end
%     
% %     lac_ci = strcat(gsm_data{ii}{LAC}, '-', ...
% %                 gsm_data{ii}{CI});
% %             
% %     latitude = gsm_data{ii}{LATITUDE};
% %     longitude = gsm_data{ii}{LONGITUDE};
%     
% %     if isempty(latitude) || isempty(longitude) ...
% %             || length(latitude)<=4 || length(longitude)<=4
% %         continue;
% %     end
% 
% %     latitude = dms2deg(gsm_data{ii}{LATITUDE});
% %     longitude = dms2deg(gsm_data{ii}{LONGITUDE});
% %     
%     
%     
%     index = index+1;
%     cell_data{index, 1} = lac_ci;
%     cell_data{index, 2} = [latitude, longitude];
% %     cell_data{ii, 3} = longitude;
% end

count = 0;
cell_data = wx_cell_info;
for jj=1:size(t_real_radius, 1)
    if t_real_radius(jj)~=0
        count = count+1;
        if isempty(cell_data{jj, 3})
            disp(['invalid data in index: ', num2str(jj)]);
            break;
        end
    else
%         cell_data{jj, 3} = [];
        if ~isempty(cell_data{jj, 3})
            disp(['invalid data in index : ', num2str(jj)]);
            break;
        end
    end
end

% t_real_radius_new = t_real_radius;
% r_real_radius_new = r_real_radius;
% 
% index = 0;
% cell_data = rogers_cell_info;
% for jj=1:size(t_real_radius, 1)
%     if t_real_radius(jj)~=0
%         t_real_radius_new(jj) = cell_data{jj, 4};
%         index = index+1;
%         r_real_radius_new(index) = cell_data{jj, 4};
%         
%         rogers_cell{jj, 11} = cell_data{jj, 4};
%     end
% end