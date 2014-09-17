% input :
%   cell_enu_angle - all the cell enu with angle info
%                  - [enu_lat, enu_long, start_angle, stop_angle, bspwr, accmin]
%   neighbour_count - how many neighbour cells be accounted
%   isd_clearance - the very closest neighbour be removed
%   toh_model - the trained OH model [A, B]
% output:
%   [radiusMax, radiusMean, radiusMin, radiusAngle]
function radius = isdAlgorithm(cell_enu_angle, neighbour_count, isd_clearance, toh_model)

idx_lat = 1;
idx_long = 2;
idx_start_angle = 3;
idx_stop_angle = 4;
idx_bspwr = 5;
idx_accmin = 6;
idx_dist = 7;

false = 0;
true = 1;

min_radius = 100;

radius = [];

for ii = 1 : length(cell_enu_angle)
    cell_enu_angle_dist = [];
    
    for jj = 1 : length(cell_enu_angle)
        if ii~=jj
            distance = sqrt((cell_enu_angle(ii, idx_lat)-cell_enu_angle(jj, idx_lat)).^2 ...
                + (cell_enu_angle(ii, idx_long)-cell_enu_angle(jj, idx_long)).^2);
        else
            distance = 0;
        end
        
        cell_enu_angle_dist = [cell_enu_angle_dist; [cell_enu_angle(jj, :), distance]];
    end
    
    [~, index] = sort(cell_enu_angle_dist(:, idx_dist));
    
    cell_enu_angle_dist = cell_enu_angle_dist(index, :);
    
    % get rif off the distance which is less than isd_clearance
    index_small = 0;
    for jj = 1 : length(cell_enu_angle_dist)
        if cell_enu_angle_dist(jj, idx_dist)>isd_clearance
            index_small = jj;
            break;
        end
    end
    
    if index_small>0
        count = length(cell_enu_angle_dist)-index_small+1;
        if count>neighbour_count
            count = neighbour_count;
        end
        
        radiusMax = max(cell_enu_angle_dist(index_small:index_small+count-1, idx_dist));
        radiusMean = mean(cell_enu_angle_dist(index_small:index_small+count-1, idx_dist));
        radiusMin = min(cell_enu_angle_dist(index_small:index_small+count-1, idx_dist));
        
        found = false;
        for jj=index_small:length(cell_enu_angle_dist)
%             inArcDirection(c, cell_angle, p)
              center.x = cell_enu_angle(ii, idx_lat);
              center.y = cell_enu_angle(ii, idx_long);
              
              p.x = cell_enu_angle_dist(jj, idx_lat);
              p.y = cell_enu_angle_dist(jj, idx_long);
              
              cell_angle = [cell_enu_angle(ii, idx_start_angle), ...
                  cell_enu_angle(ii, idx_stop_angle)];
              
              if inArcDirection(center, cell_angle, p)
                  % the cell is in direction of the cur cell
                  found = true;
                  
                  % change the cell angle for the cell
                  cell_angle = [cell_enu_angle_dist(jj, idx_start_angle), ...
                      cell_enu_angle_dist(jj, idx_stop_angle)];
                  
                  if inArcDirection(p, cell_angle, center)
                    % the cur cell is also in direction of the cell
                    % so, the radius should be ISD/2
                    radiusAngle = cell_enu_angle_dist(jj, idx_dist)/2;
                  else
                    radiusAngle = cell_enu_angle_dist(jj, idx_dist);
                  end
                  
                  break;
              end
        end
            
        if ~found
            % radiusAngle OH model 
            radiusAngle = tohFunc(cell_enu_angle(ii, :), toh_model);
        end
    elseif length(cell_enu_angle_dist)>2 %very closed
        [radiusMax, radiusMean, radiusMin, radiusAngle] = deal(min_radius);
    else
        % OH model
        radiusAngle = tohFunc(cell_enu_angle(ii, :), toh_model);
        
        [radiusMax, radiusMean, radiusMin] = deal(radiusAngle);
    end     
    
    radius = [radius; radiusMax, radiusMean, radiusMin, radiusAngle];
end

end


function radius = tohFunc(cell_enu_angle, toh_model)

max_radius = 30000;

%[A, B, radius] = tohPropagationAlgorithm(cell_radius_power, toh_model)
if ~isempty(toh_model)
    cell_radius_power = [cell_enu_angle(:, idx_bsspwr), ...
        cell_enu_angle(:, idx_accmin), 0];
    [~, ~, radius] = tohPropagationAlgorithm(cell_radius_power, toh_model);
else
    radius = max_radius;
end

radius = min(max_radius, radius);

end