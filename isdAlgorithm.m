% input :
%   cell_enu_angle - all the cell enu with angle info
%                  - [enu_lat, enu_long, start_angle, stop_angle, bspwr, accmin, real_radius]
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

isd_radius = zeros(length(cell_enu_angle));

% compute all the isd distances, and store in isd_radius
load rogers_isd_radius;
isd_radius = rogers_isd_radius_all;

% for ii = 1 : length(cell_enu_angle)
% %     cell_enu_angle_dist = [];
%     
%     for jj = 1 : length(cell_enu_angle)
%         if ii~=jj
%             if isd_radius(jj,ii)~=0
%                 isd_radius(ii, jj)= isd_radius(jj, ii);
%             else
%                 isd_radius(ii, jj) = sqrt((cell_enu_angle(ii, idx_lat)-cell_enu_angle(jj, idx_lat)).^2 ...
%                     + (cell_enu_angle(ii, idx_long)-cell_enu_angle(jj, idx_long)).^2);
%             end
%         end
%     end
% end
   

for ii=1:length(cell_enu_angle)
    
    if cell_enu_angle(ii, idx_dist)~=0 %should be compute the radius
        [~, index] = sort(isd_radius(:, ii));
    
        cell_enu_angle_dist = cell_enu_angle(index, :);
        cell_enu_angle_dist(:, idx_dist) = isd_radius(index, ii);
    
        % get rif off the distance which is less than isd_clearance
        index_small = 0;
        index_max = 0;
        unq_radius = [];
        for jj = 1 : length(cell_enu_angle_dist)
            if cell_enu_angle_dist(jj, idx_dist)>isd_clearance
                if ~index_small
                    index_small = jj;
                end
                
                index_max = jj;
                
                unq_radius = unique(cell_enu_angle_dist(index_small:index_max, idx_dist));
                if size(unq_radius, 1)>=neighbour_count
                    break;
                end
            end
        end
    
        if index_small>0
%             count = length(cell_enu_angle_dist)-index_small+1;
%             if count>neighbour_count
%                 count = neighbour_count;
%             end
            
            radiusMax = max(unq_radius);
            radiusMean = mean(unq_radius);
            radiusMin = min(unq_radius);
            
            found = false;
            % for jj=index_small:length(cell_enu_angle_dist)
            for jj=index_small:index_max
                
                %inArcDirection(c, cell_angle, p)
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

end


function radius = tohFunc(cell_enu_angle, toh_model)

idx_bspwr = 5;
idx_accmin = 6;

max_radius = 30000;

%[A, B, radius] = tohPropagationAlgorithm(cell_radius_power, toh_model)
if ~isempty(toh_model)
    cell_radius_power = [cell_enu_angle(:, idx_bspwr), ...
        cell_enu_angle(:, idx_accmin), 0];
    [~, ~, radius] = tohPropagationAlgorithm(cell_radius_power, toh_model);
else
    radius = max_radius;
end

radius = min(max_radius, radius);

end