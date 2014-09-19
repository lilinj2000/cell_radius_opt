
% close all;
% clear all;
% clc;
% 
% load rogers_radius;
% load rogers_voronoi;
% load rogers_isd_radius;

show_isd = true;
show_voronoi = false;

for index = 49 : 50
    cell_enu = r_cell_enu(index, :);
    cell_angle = r_cell_angle(index, :);
    
    if show_voronoi
        
        voro_smax = 5;
        radius = rogers_voro_radius(index, voro_smax);
    elseif show_isd
        isd_max = 1;
        radius = rogers_isd_radius(index, isd_max);
    end
        
    d_cell_enu_angle_radius = [cell_enu, cell_angle, radius];

    t_cell_enu_angle = [t_cell_enu, t_cell_angle];
        
        
    if show_voronoi
        [segs, vertexs, nbr_points] = rogers_voro.findResultInfo(cell_enu);
        
        drawCellVoroGeo(d_cell_enu_angle_radius, vertexs, nbr_points, segs, t_cell_enu_angle);

    end
    
    if show_isd
       
        p_loc = find(ismember(t_cell_enu, cell_enu, 'rows')>0);
%         if length(p_loc)>1 % multi cell, then compare an
%             p_cell_angle = t_cell_angle(p_loc, :);
%             p_loc = find(ismember(p_cell_angle, cell_angle, 'rows')>0);
%         end
        
        
        isd_radius = rogers_isd_radius_all(p_loc(1), :);
        drawCellISDGeo(d_cell_enu_angle_radius, t_cell_enu_angle, isd_radius, 5, 50);
    end
end