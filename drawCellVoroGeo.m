function drawCellVoroGeo(cell_enu, cell_angle, cell_radius, vertexs, nbr_points, segs, t_cell_enu, t_cell_angle)

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

figure;
hold on;

% draw point
plot(p.x, p.y, '*');

% draw vertex
plot(vertexs(:, 1), vertexs(:, 2), '^');

% draw neighbour site points
plot(nbr_points(:, 1), nbr_points(:, 2), '*');

% draw segs
for ii=1:size(segs, 2)
    drawLine(segs(ii).start_p, segs(ii).end_p);
end


% draw arc
[start_angle, stop_angle] = changeGeoAngle(start_angle, stop_angle);
drawArc(p, radius, start_angle, stop_angle, 1);

% draw circle
[X, Y ] = drawCircle(p.x, p.y, radius);
plot(X, Y, '--');

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
                    nbr_start_angle = t_cell_enu_angle(jj, idx_start_angle);
                    nbr_stop_angle = t_cell_enu_angle(jj, idx_stop_angle);
                end
                [nbr_start_angle, nbr_stop_angle] = changeGeoAngle(nbr_start_angle, ...
                    nbr_stop_angle);
                
                drawArc(nbr_point, radius/4, nbr_start_angle, nbr_stop_angle, 1);
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