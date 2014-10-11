% input :
%   cell_enu    - [enu_lat, enu_long]
%   cell_angle  - [start_angle, stop_angle]
%   cell_power  - [bspwr, accmin] 
%   cell_real_radius - [real_radius]
%   neighbour_count - how many neighbour cells be accounted
%   isd_clearance - the very closest neighbour be removed
%   toh_model - the trained OH model [A, B]
% output:
%   [radiusMax, radiusMean, radiusMin, radiusAngle]
function radius = isdAlgorithm(cell_enu, cell_angle, cell_power, cell_real_radius, neighbour_count, isd_clearance, toh_model)

idx_lat = 1;
idx_long = 2;

idx_start_angle = 1;
idx_stop_angle = 2;
idx_bspwr = 1;
idx_accmin = 2;

idx_dist = 1;

log = log4m.getLogger('isd.log');
log.setCommandWindowLevel(log.OFF);

min_radius = 100;

radius = [];

isd_radius = zeros(size(cell_enu, 1));

% compute all the isd distances, and store in isd_radius
% load turkey_isd_radius;
% isd_radius = turkey_isd_radius_all;

load rogers_isd_radius;
isd_radius = rogers_isd_radius_all;

% for ii = 1 : size(cell_enu, 1)
% %     cell_enu_angle_dist = [];
%     
%     log.info('isdAlgorithm()', ['compute no ', num2str(ii)]);
% 
%     for jj = 1 : size(cell_enu, 1)
%         if ii~=jj
%             if isd_radius(jj,ii)~=0
%                 isd_radius(ii, jj)= isd_radius(jj, ii);
%             else
%                 isd_radius(ii, jj) = sqrt((cell_enu(ii, idx_lat)-cell_enu(jj, idx_lat)).^2 ...
%                     + (cell_enu(ii, idx_long)-cell_enu(jj, idx_long)).^2);
%             end
%         end
%     end
% end
   

for ii=1:size(cell_enu,1)
    
    if cell_real_radius(ii, idx_dist)~=0 %should be compute the radius
        
        log.info('isdAlgorithm()', ['compute radius ', num2str(ii)]);
        
        [~, index] = sort(isd_radius(:, ii));
    
        cell_enu_tmp = cell_enu(index, :);
        
        if ~isempty(cell_angle)
            cell_angle_tmp = cell_angle(index, :);
        end
        
%         if ~isempty(cell_power)
%             cell_power_tmp = cell_power(index, :);
%         end
        
        cell_dist = isd_radius(index, ii);
    
        % get rif off the distance which is less than isd_clearance
        index_small = 0;
        index_max = 0;
        unq_radius = [];
        for jj = 1 : size(cell_dist, 1)
            if cell_dist(jj, idx_dist)>isd_clearance
                if ~index_small
                    index_small = jj;
                end
                
                index_max = jj;
                
                unq_radius = unique(cell_dist(index_small:index_max, idx_dist));
                if size(unq_radius, 1)>=neighbour_count
                    break;
                end
            end
        end
    
        if index_small>0

            radiusMax = max(unq_radius);
            radiusMean = mean(unq_radius);
            radiusMin = min(unq_radius);
            
            if isempty(cell_angle)
                radiusAngle = radiusMin;
            else
                found = false;
                % for jj=index_small:length(cell_enu_angle_dist)
                for jj=index_small:index_max
                    
                    %inArcDirection(c, cell_angle, p)
                    center.x = cell_enu(ii, idx_lat);
                    center.y = cell_enu(ii, idx_long);
                    
                    if found %% check the same point
                        if ~(p.x==cell_enu_tmp(jj, idx_lat) ...
                            && p.y==cell_enu_tmp(jj, idx_long))
                            break;
                        end
                    end
                    
                    p.x = cell_enu_tmp(jj, idx_lat);
                    p.y = cell_enu_tmp(jj, idx_long);
                    
                    cur_cell_angle = [cell_angle(ii, idx_start_angle), ...
                        cell_angle(ii, idx_stop_angle)];
                    
                    
                    if inArcDirection(center, cur_cell_angle, p)
                        % the cell is in direction of the cur cell
                        found = true;
                        
                        % change the cell angle for the cell
                        cur_cell_angle = [cell_angle_tmp(jj, idx_start_angle), ...
                            cell_angle_tmp(jj, idx_stop_angle)];
                        
                        if inArcDirection(p, cur_cell_angle, center)
                            % the cur cell is also in direction of the cell
                            % so, the radius should be ISD/2
                            radiusAngle = cell_dist(jj, idx_dist)/2;
                            break;
                        else
                            radiusAngle = cell_dist(jj, idx_dist);
                            continue; % check the next point whether the same point
                        end
                    end
                end
                
                if ~found
                    % radiusAngle OH model
                    radiusAngle = tohFunc(cell_power(ii, :), toh_model);
                end
            end
        elseif length(cell_dist)>2 %very closed
            [radiusMax, radiusMean, radiusMin, radiusAngle] = deal(min_radius);
        else
            % OH model
            radiusAngle = tohFunc(cell_power(ii, :), toh_model);
            
            [radiusMax, radiusMean, radiusMin] = deal(radiusAngle);
        end
    
        radius = [radius; radiusMax, radiusMean, radiusMin, radiusAngle];
    end
end

end


function radius = tohFunc(cell_power, toh_model)

max_radius = 30000;

if isempty(cell_power)
    radius = max_radius;
    return ;
end
                        
idx_bspwr = 1;
idx_accmin = 2;


%[A, B, radius] = tohPropagationAlgorithm(cell_radius_power, toh_model)
if ~isempty(toh_model)
    cell_radius_power = [cell_power(:, idx_bspwr), ...
        cell_power(:, idx_accmin), 0];
    [~, ~, radius] = tohPropagationAlgorithm(cell_radius_power, toh_model);
else
    radius = max_radius;
end

radius = min(max_radius, radius);

end