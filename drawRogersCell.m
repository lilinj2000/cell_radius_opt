
% close all;
% clear all;
% clc;
% 


% load india_radius;
% load india_voronoi;
% load india_isd_radius;
% load wx_radius;
% load wx_voronoi;
% load wx_isd_radius;


show_isd = true;
show_voronoi = true;

by_cellid =  false;

%% rogers
% load rogers_radius;
% load rogers_voronoi;
% load rogers_isd_radius;
% cell_info_all = rogers_cell_info;
% 
% voro_radius = rogers_voro_radius;
% isd_radius = rogers_isd_radius;
% isd_radius_all = rogers_isd_radius_all;
% v = rogers_voro;



%% turkey
% load turkey_radius;
% load turkey_voronoi;
% load turkey_isd_radius;
% cell_info_all = turkey_cell_info;
% 
% voro_radius = turkey_voro_radius;
% isd_radius = turkey_isd_radius;
% isd_radius_all = turkey_isd_radius_all;
% v = turkey_voro;
% 
% lac_ci = '35506-58331';

%% china wuxi
% cell_info_all = wx_cell_info;
% 
% voro_radius = wx_voro_radius;
% isd_radius = wx_isd_radius;
% isd_radius_all = wx_isd_radius_all;
% v = wx_voro;
% lac_ci = '0-112';
% lac_ci = '0-113';

%% china beijing
% load beijing_radius;
% load beijing_voro;
% load beijing_isd_radius;

cell_info_all = beijing_cell_total_new;
% voro_radius = beijing_voro_radius;
% isd_radius = beijing_isd_radius;
voro_radius = beijing_voro_best_radius;
isd_radius = beijing_isd_best_radius;
isd_radius_all = beijing_isd_radius_all_new;
v = beijing_voro_new;
by_cellid =  true;
% lac_ci = '10687-54851';
% lac_ci = '10712-54826';

close all;

idx_real_radius = 2; % 0.95 mesrs

% 
for index=10:11
    if by_cellid
        t_lac_ci = t_cellid;
        lac_ci = r_cellid(index);
        lac_ci_str = num2str(lac_ci);
    else
        lac_ci = r_lac_ci{index, 1};
        lac_ci_str = lac_ci;
    end
    
%     index = find(ismember(r_lac_ci, lac_ci)>0);

    
    if index>0
        
        if ~isempty(r_cell_angle)
            cell_angle = r_cell_angle(index, :);
        else
            cell_angle = [];
        end
        
        if by_cellid
            
            idx_cell_start = find(cell_info_all(:, 3)==lac_ci, 1, 'first');
            idx_cell_end = find(cell_info_all(:, 3)==lac_ci, 1, 'last');
            
            
            for kk=idx_cell_start:idx_cell_end
                if cell_info_all(kk, 5)~=0
                    cell_info{1, 1} = lac_ci_str;
                    cell_info{1, 2} = cell_info_all(kk, 1:2);
                    cell_info{1, 3} = cell_info_all(kk, 5:7);
                end
            end
            
            cell_lat_long = cell_info{1, 2};
            
            idx_enu = find(ismember(unq_lat_long, cell_lat_long, 'rows')>0);
            
            cell_enu = unq_enu(idx_enu, :);
            t_cell_enu = unq_enu;
%             cell_count = idx_cell_end - idx_cell_start + 1;
%             
%             cell_info = cell(cell_count, 3);
%             cell_info(:, 1) = {lac_ci};
%             cell_info(:, 2) = num2cell(cell_info_all(idx_cell_start:idx_cell_end, 2:3), 2);
%             cell_info(:, 3) = num2cell(cell_info_all(idx_cell_start:idx_cell_end, 5:7), 2);
        else
            index_loc = find(ismember(t_lac_ci, lac_ci)>0);
            cell_info = cell_info_all(index_loc, :);
            
            cell_enu = r_cell_enu(index, :);
        end
        
        drawCellMsrGeo(cell_info, cell_angle);
        cell_radius = cell_info{1, 3}(:, idx_real_radius);
        
        if show_voronoi
            v_radius = voro_radius(index, :);
        end
        
        if show_isd
            i_radius = isd_radius(index, :);
        end
        
        %     d_cell_enu_angle_radius = [cell_enu, cell_angle, radius];
        %
        %     t_cell_enu_angle = [t_cell_enu, t_cell_angle];
        
        
        if show_voronoi
            [segs, vertexs, nbr_points] = v.findResultInfo(cell_enu);
            
            drawCellVoroGeo(lac_ci_str, cell_enu, cell_angle, cell_radius, v_radius, vertexs, nbr_points, segs, t_cell_enu, t_cell_angle);
            
        end
        
        if show_isd
            
            p_loc = find(ismember(t_cell_enu, cell_enu, 'rows')>0);
            %         if length(p_loc)>1 % multi cell, then compare an
            %             p_cell_angle = t_cell_angle(p_loc, :);
            %             p_loc = find(ismember(p_cell_angle, cell_angle, 'rows')>0);
            %         end
            
            
            t_isd_radius = isd_radius_all(p_loc(1), :);
            drawCellISDGeo(lac_ci_str, cell_enu, cell_angle, cell_radius, i_radius, t_cell_enu, t_cell_angle, t_isd_radius, 5, 50);
        end
    end
end