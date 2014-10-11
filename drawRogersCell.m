
% close all;
% clear all;
% clc;
% 
% load rogers_radius;
% load rogers_voronoi;
% load rogers_isd_radius;
% load turkey_radius;
% load turkey_voronoi;
% load turkey_isd_radius;
% load india_radius;
% load india_voronoi;
% load india_isd_radius;
% load wx_radius;
% load wx_voronoi;
% load wx_isd_radius;


show_isd = true;
show_voronoi = true;

%% rogers
cell_info_all = rogers_cell_info;

voro_radius = rogers_voro_radius;
isd_radius = rogers_isd_radius;
isd_radius_all = rogers_isd_radius_all;
v = rogers_voro;

% lac_ci = '4300-25077';
% lac_ci = '17500-33411';
% lac_ci = '17500-33418';
% lac_ci = '17500-33419';
% lac_ci = '4400-28201';
% lac_ci = '4400-37528';

%% turkey
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
% % lac_ci = '0-112';
% lac_ci = '0-113';

%% china beijing
% cell_info_all = beijing_cell_info;
% voro_radius = beijing_voro_radius;
% isd_radius = beijing_isd_radius;
% isd_radius_all = beijing_isd_radius_all;
% v = beijing_voro;
% lac_ci = '10687-54851';
% lac_ci = '10712-54826';

close all;

% 
for index=255:256
    lac_ci = r_lac_ci{index, 1};
    
%     index = find(ismember(r_lac_ci, lac_ci)>0);

    
    if index>0
        
        if ~isempty(r_cell_angle)
            cell_angle = r_cell_angle(index, :);
        else
            cell_angle = [];
        end
        
        index_loc = find(ismember(t_lac_ci, lac_ci)>0);
        cell_info = cell_info_all(index_loc, :);
        drawCellMsrGeo(cell_info, cell_angle);
        
        cell_enu = r_cell_enu(index, :);
        
        
        
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
            
            drawCellVoroGeo(lac_ci, cell_enu, cell_angle, v_radius, vertexs, nbr_points, segs, t_cell_enu, t_cell_angle);
            
        end
        
        if show_isd
            
            p_loc = find(ismember(t_cell_enu, cell_enu, 'rows')>0);
            %         if length(p_loc)>1 % multi cell, then compare an
            %             p_cell_angle = t_cell_angle(p_loc, :);
            %             p_loc = find(ismember(p_cell_angle, cell_angle, 'rows')>0);
            %         end
            
            
            t_isd_radius = isd_radius_all(p_loc(1), :);
            drawCellISDGeo(lac_ci, cell_enu, cell_angle, i_radius, t_cell_enu, t_cell_angle, t_isd_radius, 5, 50);
        end
    end
end