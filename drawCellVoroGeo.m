function drawCellVoroGeo(lac_ci, cell_enu, cell_angle, cell_radius, v_radius, vertexs, nbr_points, segs, t_cell_enu, t_cell_angle)

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

idx_v_max = 1;
idx_v_mean = 2;
idx_v_min = 3;
idx_v_angle = 4;
idx_s_max = 5;
idx_s_mean = 6;
idx_s_min = 7;
idx_s_angle = 8;

radius = v_radius(idx_s_max);
[start_angle, stop_angle] = changeGeoAngle(start_angle, stop_angle);

show_v = true;
show_v_angle = false;
show_s = false;
show_s_angle = false;

steps = [];
if show_v
    steps = [steps; 1]; 
end

if show_v_angle
    steps = [steps; 2];
end

if show_s
    steps = [steps; 3];
end

if show_s_angle
    steps = [steps; 4];
end

for step=steps
    figure;
    hold on;
    
    % draw point
    plot(p.x, p.y, '*');
    text(p.x+5, p.y+5, lac_ci);
    
    % draw vertex
    plot(vertexs(:, 1), vertexs(:, 2), '^');
    
    % draw neighbour site points
    plot(nbr_points(:, 1), nbr_points(:, 2), '*');
    
    % draw segs
    for ii=1:size(segs, 2)
        drawLine(segs(ii).start_p, segs(ii).end_p);
    end


    % draw arc
    drawArc(p, radius, start_angle, stop_angle, 1);

    % draw circle
    
    [X, Y ] = drawCircle(p.x, p.y, cell_radius);
    plot(X, Y, 'k--', 'LineWidth',2);
        
    if step==1
        [X, Y ] = drawCircle(p.x, p.y, v_radius(idx_v_max));
        plot(X, Y, '--');
        
        [X, Y ] = drawCircle(p.x, p.y, v_radius(idx_v_mean));
        plot(X, Y, 'g--');
        
        [X, Y ] = drawCircle(p.x, p.y, v_radius(idx_v_min));
        plot(X, Y, 'r--');
        
        
        
    elseif step==2
        
        [X, Y ] = drawCircle(p.x, p.y, v_radius(idx_v_angle));
        plot(X, Y, 'c--');
        
        
    elseif step==3
        [X, Y ] = drawCircle(p.x, p.y, v_radius(idx_s_max));
        plot(X, Y, '--');
        
        [X, Y ] = drawCircle(p.x, p.y, v_radius(idx_s_mean));
        plot(X, Y, 'g--');
        
        [X, Y ] = drawCircle(p.x, p.y, v_radius(idx_s_min));
        plot(X, Y, 'r--');
        
    elseif step==4
        [X, Y ] = drawCircle(p.x, p.y, v_radius(idx_s_angle));
        plot(X, Y, 'c--');
    end

    % draw arc for the same direction neighbour points
    for ii=1:size(nbr_points, 1)
        
        nbr_point.x = nbr_points(ii, 1);
        nbr_point.y = nbr_points(ii, 2);
        
        
        if inArcDirection(p, cell_angle(1, idx_start_angle:idx_stop_angle), nbr_point)
            for jj=1:size(t_cell_enu, 1)
                if t_cell_enu(jj, idx_lat)==nbr_point.x ...
                        && t_cell_enu(jj, idx_long)==nbr_point.y
                    if isempty(t_cell_angle)
                        continue;
                    else
                        nbr_start_angle = t_cell_angle(jj, idx_start_angle);
                        nbr_stop_angle = t_cell_angle(jj, idx_stop_angle);
                    end
                    [nbr_start_angle, nbr_stop_angle] = changeGeoAngle(nbr_start_angle, ...
                        nbr_stop_angle);
                    
                    distance = sqrt((nbr_point.x-p.x).^2 + (nbr_point.y-p.y).^2);
                    drawArc(nbr_point, distance/2, nbr_start_angle, nbr_stop_angle, 1);
                end
            end
        end
    end


% [x_min, x_max] = deal(min(nbr_points(:, 1)), max(nbr_points(:, 1)));
% [y_min, y_max] = deal(min(nbr_points(:, 2)), max(nbr_points(:, 2)));
% 
% deltaX = (x_max - x_min)*0.25 ;
% deltaY = (y_max - y_min)*0.25 ;


% axis([x_min-deltaX, x_max+deltaX, y_min-deltaY, y_max+deltaY]);

    axis equal;
end

end