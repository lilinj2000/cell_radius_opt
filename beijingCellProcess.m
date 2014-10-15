%% beijing cell process

% idx_lat = 1;
% idx_long = 2;
% idx_cellid = 14;
% cell_total = beijing_cell_total(:, [idx_lat, idx_long, idx_cellid]);
% 
% cell_total = unique(cell_total, 'rows');
% 
% cell_total = sortrows(cell_total, [3, 1, 2]);
% 
% 
% cell_total(:, 4) = 1:size(cell_total, 1);
% 
% for ii=1:size(cell_total, 1)-1
%     
%     if cell_total(ii, 3)==cell_total(ii+1, 3)
%         cell_total(ii+1, 4) = cell_total(ii, 4);
%     end
%     
% end
% 
% 
% bj_msr = measurement(:, [2, 18, 17]);
% 
% bj_msr = sortrows(bj_msr, [1, 2, 3]);

% cell_total(:, [5, 6, 7]) = 0;
% 
% unq_cellid = sort(unique(cell_total(:, 3)));
% 
% for ii=1:size(unq_cellid, 1)
%     cellid = unq_cellid(ii);
%     
%     idx_cell_start = find(cell_total(:, 3)==cellid, 1, 'first');
%     idx_cell_end = find(cell_total(:, 3)==cellid, 1, 'last');
%     cell_count = idx_cell_end - idx_cell_start + 1;
%     
%     idx_msr_start = find(bj_msr(:, 1)==cellid, 1, 'first'); 
%     idx_msr_end = find(bj_msr(:, 1)==cellid, 1, 'last');
%     msr_count = idx_msr_end - idx_msr_start + 1;
%     
%     
%     %% compute the distance between msr and cell
%     dist = zeros(msr_count, cell_count);
%     for jj=1:cell_count
%         ref_lat_long = [cell_total(idx_cell_start+jj-1, 1), cell_total(idx_cell_start+jj-1, 2)];
%         ref_lat_long = deg2rad(ref_lat_long);
%         
%         msr_lat_long = deg2rad(bj_msr(idx_msr_start:idx_msr_end, 2:3));
%         
%         msr_enu = convertlatlong2enu(msr_lat_long, ref_lat_long);
%         
%         for kk=1:msr_count
%             dist(kk, jj) = sqrt(msr_enu(kk, 1).^2+msr_enu(kk, 2).^2);
%         end
%     end
%     
%     % set the non min dist to zero
%     for jj=1:size(dist, 1)
%         min_dist = min(dist(jj, :));
%         
%         idx = find(dist(jj, :)~=min_dist>0);
%         dist(jj, idx) = 0;
%     end
%     
%     for jj=1:cell_count
%         cell_radius = dist(:, jj);
%         
%         cell_radius(find(cell_radius==0)) = [];
%         
%         if ~isempty(cell_radius) 
%             cell_radius = sort(cell_radius);
%             cell_total(idx_cell_start+jj-1, 5) = cell_radius(round(size(cell_radius, 1)*0.8));
%             cell_total(idx_cell_start+jj-1, 6) = cell_radius(round(size(cell_radius, 1)*0.95));
%             cell_total(idx_cell_start+jj-1, 7) = cell_radius(size(cell_radius, 1));
%         end
%     end
% end


% cell_total = beijing_cell_total_new ;
% 
% idx_lat_long = [1, 2];
% idx_cellid = 3;
% 
% idx_radius = [5, 6, 7];
% 
% t_lat_long = zeros(size(cell_total, 1), 2);
% r_lat_long = [];
% t_cellid = zeros(size(cell_total, 1), 1);
% r_cellid = [];
% t_real_radius = zeros(size(cell_total, 1), 3);
% r_real_radius = [];
% 
% for ii=1:size(cell_total, 1)
%     t_lat_long(ii, :) = cell_total(ii, idx_lat_long);
%     t_cellid(ii, :) = cell_total(ii, idx_cellid);
%     t_real_radius(ii, :) = cell_total(ii, idx_radius);
%     
%     if ~isempty(find(cell_total(ii, idx_radius)~=0>0))
%         r_lat_long = [r_lat_long; t_lat_long(ii, :)];
%         r_cellid = [r_cellid; t_cellid(ii, :)];
%         r_real_radius = [r_real_radius; t_real_radius(ii, :)];
%     end
% end


%% compute isd radius
% min_radius = 100;
% max_radius = 30000;
% beijing_isd_radius = zeros(size(r_real_radius, 1), 3);
% 
% for ii=1:size(r_real_radius, 1)
%     lat_long = r_lat_long(ii, :);
%     index = find(ismember(unq_lat_long, lat_long, 'rows')>0);
%     
%     cell_radius = beijing_isd_radius_all_new(:, index);
%     cell_radius = sort(cell_radius);
%     
%     idx_start = find(cell_radius>50, 1, 'first');
%     unq_radius = [];
%     for jj=idx_start:size(cell_radius, 1)
%         unq_radius = unique(cell_radius(idx_start:jj));
%         
%         if size(unq_radius,1)>=5
%             break;
%         end
%     end
%     
%     if isempty(unq_radius)
%         [radius_max, radius_mean, radius_min] = deal(min_radius);
%     else
%         radius_max = max(unq_radius);
%         radius_mean = mean(unq_radius);
%         radius_min = min(unq_radius);
%     end
%     
%     beijing_isd_radius(ii, :) = [radius_max, radius_mean, radius_min];
% end

%% compute voro radius
cell_angle = [];
cell_power = [];
isd_clearance = 50;
v = beijing_voro_new;
beijing_voro_radius = zeros(size(r_real_radius, 1), 8);
min_radius = 100;
max_radius = 30000;
for ii=1:size(r_real_radius, 1)
    
    lat_long = r_lat_long(ii, :);
    index = find(ismember(unq_lat_long, lat_long, 'rows')>0);
    
    enu_latlong = unq_enu(index, :);
    p.x = enu_latlong(1, 1);
    p.y = enu_latlong(1, 2);
    
    vradius = v.computeVRadius(p, cell_angle, cell_power, isd_clearance, max_radius, min_radius, toh_model);
    
    sradius = v.computeSRadius(p, cell_angle, cell_power, unq_enu, cell_angle, isd_clearance, max_radius, min_radius, toh_model);
    
    beijing_voro_radius(ii, :) = [vradius, sradius];
end

