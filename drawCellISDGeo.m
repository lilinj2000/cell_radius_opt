function drawCellISDGeo(cell_enu, cell_angle, cell_radius, t_cell_enu, t_cell_angle, isd_radius, nbr_count, isd_clearance)

idx_lat = 1;
idx_long = 2;
idx_start_angle = 1;
idx_stop_angle = 2;
idx_radius = 1;

p.x = cell_enu(idx_lat);
p.y = cell_enu(idx_long);

if isempty(cell_angle)
    cell_angle = [0, 360];
end

start_angle = cell_angle(idx_start_angle);
stop_angle = cell_angle(idx_stop_angle);


radius = cell_radius(idx_radius);

[~, index] = sort(isd_radius(:));
t_cell_enu = t_cell_enu(index, :);
if ~isempty(t_cell_angle)
    t_cell_angle = t_cell_angle(index, :);
end
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
    if cur_cell>size(t_cell_enu, 1)
        break;
    end
    
    nbr_point.x = t_cell_enu(cur_cell, idx_lat);
    nbr_point.y = t_cell_enu(cur_cell, idx_long);
    
   
    if isd_radius(cur_cell)>isd_clearance

         plot(nbr_point.x, nbr_point.y, '*');
    
         if inArcDirection(p, cell_angle(1, idx_start_angle:idx_stop_angle), nbr_point)
             for jj=1:size(t_cell_enu, 1)
                 if t_cell_enu(jj, idx_lat)==nbr_point.x ...
                         && t_cell_enu(jj, idx_long)==nbr_point.y
                     if ~isempty(t_cell_angle)
                        nbr_start_angle = t_cell_enu_angle(jj, idx_start_angle);
                        nbr_stop_angle = t_cell_enu_angle(jj, idx_stop_angle);
                     else
                         continue;
                     end
                     
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