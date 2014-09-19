function drawCellISDGeo(cell_enu_angle_radius, t_cell_enu_angle, isd_radius, nbr_count, isd_clearance)

idx_lat = 1;
idx_long = 2;
idx_start_angle = 3;
idx_stop_angle = 4;
idx_radius = 5;

p.x = cell_enu_angle_radius(idx_lat);
p.y = cell_enu_angle_radius(idx_long);
start_angle = cell_enu_angle_radius(idx_start_angle);
stop_angle = cell_enu_angle_radius(idx_stop_angle);
radius = cell_enu_angle_radius(idx_radius);

[~, index] = sort(isd_radius(:));
t_cell_enu_angle = t_cell_enu_angle(index, :);
isd_radius = isd_radius(index);

figure;
hold on;

% draw point
plot(p.x, p.y, '*');

% draw arc
[start_angle, stop_angle] = changeGeoAngle(start_angle, stop_angle);
drawArc(p, radius, start_angle, stop_angle, 1);

% draw circle
[X, Y ] = drawCircle(p.x, p.y, radius);
plot(X, Y, '--');

% draw nbr points
cur_cell = 0;
radius_set = [];
while true
    cur_cell = cur_cell + 1;
    if cur_cell>size(t_cell_enu_angle, 1)
        break;
    end
    
    nbr_point.x = t_cell_enu_angle(cur_cell, idx_lat);
    nbr_point.y = t_cell_enu_angle(cur_cell, idx_long);
    
   
    if isd_radius(cur_cell)>isd_clearance

         plot(nbr_point.x, nbr_point.y, '*');
    
         if inArcDirection(p, cell_enu_angle_radius(1, idx_start_angle:idx_stop_angle), nbr_point)
             for jj=1:size(t_cell_enu_angle, 1)
                 if t_cell_enu_angle(jj, idx_lat)==nbr_point.x ...
                         && t_cell_enu_angle(jj, idx_long)==nbr_point.y
                     nbr_start_angle = t_cell_enu_angle(jj, idx_start_angle);
                     nbr_stop_angle = t_cell_enu_angle(jj, idx_stop_angle);
                     [nbr_start_angle, nbr_stop_angle] = changeGeoAngle(nbr_start_angle, ...
                         nbr_stop_angle);
                     
                     drawArc(nbr_point, radius/4, nbr_start_angle, nbr_stop_angle, 1);
                 end
             end
         end
         
         if any(radius_set==isd_radius(cur_cell))
            continue;
         end
        
        radius_set = [radius_set; isd_radius(cur_cell)];
        
        if size(radius_set, 1)>=nbr_count
            break;
        end
    else
        if nbr_point.x==p.x && nbr_point.y==p.y
            continue;
        end
        
        plot(nbr_point.x, nbr_point.y, 'o');
    end
    
end


end