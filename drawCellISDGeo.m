function drawCellISDGeo(lac_ci, cell_enu, cell_angle, cell_radius, i_radius, t_cell_enu, t_cell_angle, isd_radius, nbr_count, isd_clearance)

idx_lat = 1;
idx_long = 2;
idx_start_angle = 1;
idx_stop_angle = 2;
% idx_radius = 1;

idx_max = 1;
idx_mean = 2;
idx_min = 3;
idx_angle = 4;

p.x = cell_enu(idx_lat);
p.y = cell_enu(idx_long);

if isempty(cell_angle)
    cell_angle = [0, 360];
end

start_angle = cell_angle(idx_start_angle);
stop_angle = cell_angle(idx_stop_angle);


radius = i_radius(idx_max);

[~, index] = sort(isd_radius(:));
t_cell_enu = t_cell_enu(index, :);
if ~isempty(t_cell_angle)
    t_cell_angle = t_cell_angle(index, :);
end
isd_radius = isd_radius(index);

[start_angle, stop_angle] = changeGeoAngle(start_angle, stop_angle);


show_i = true;
show_i_angle = false;

steps = [];
if show_i
    steps = [steps; 1];
end

if show_i_angle
    steps = [steps; 2];
end

for step = steps
    figure;
    hold on;
    
    % draw point
    plot(p.x, p.y, '*');
    text(p.x+5, p.y+5, lac_ci);
    
    % draw arc
    drawArc(p, radius, start_angle, stop_angle, 1);
    
    % draw circle
    [X, Y ] = drawCircle(p.x, p.y, cell_radius);
    plot(X, Y, 'k--', 'LineWidth',2);
    
    
    if step==1
        [X, Y ] = drawCircle(p.x, p.y, i_radius(idx_max));
        plot(X, Y, '--');
        
        [X, Y ] = drawCircle(p.x, p.y, i_radius(idx_mean));
        plot(X, Y, 'g--');
        
        [X, Y ] = drawCircle(p.x, p.y, i_radius(idx_min));
        plot(X, Y, 'r--');
        
    elseif step==2
        [X, Y] = drawCircle(p.x, p.y, i_radius(idx_angle));
        plot(X, Y, 'c--');
    end
    
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
                            nbr_start_angle = t_cell_angle(jj, idx_start_angle);
                            nbr_stop_angle = t_cell_angle(jj, idx_stop_angle);
                        else
                            continue;
                        end
                        
                        [nbr_start_angle, nbr_stop_angle] = changeGeoAngle(nbr_start_angle, ...
                            nbr_stop_angle);
                        
                        drawArc(nbr_point, isd_radius(cur_cell)/2, nbr_start_angle, nbr_stop_angle, 1);
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


end