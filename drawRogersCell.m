
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
show_voronoi = false;

voro_radius = rogers_voro_radius;
isd_radius = rogers_isd_radius;
isd_radius_all = rogers_isd_radius_all;
v = rogers_voro;

for index = 35 : 35
    cell_enu = r_cell_enu(index, :);
    
    if ~isempty(r_cell_angle)
        cell_angle = r_cell_angle(index, :);
    else
        cell_angle = [];
    end
    
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
        
        drawCellVoroGeo(cell_enu, cell_angle, v_radius, vertexs, nbr_points, segs, t_cell_enu, t_cell_angle);

    end
    
    if show_isd
       
        p_loc = find(ismember(t_cell_enu, cell_enu, 'rows')>0);
%         if length(p_loc)>1 % multi cell, then compare an
%             p_cell_angle = t_cell_angle(p_loc, :);
%             p_loc = find(ismember(p_cell_angle, cell_angle, 'rows')>0);
%         end
        
        
        t_isd_radius = isd_radius_all(p_loc(1), :);
        drawCellISDGeo(cell_enu, cell_angle, i_radius, t_cell_enu, t_cell_angle, t_isd_radius, 5, 50);
    end
end