% split rogers cell

% fetch the lac ci info
% idx_lac_ci = 1;
% idx_real_dist = 11;
% t_lac_ci = cell(0);
% r_lac_ci = cell(0);
% for ii=1:size(rogers_cell)
%     t_lac_ci(size(t_lac_ci, 1)+1, 1) = cellstr(rogers_cell{ii, idx_lac_ci}); 
%     if rogers_cell{ii, idx_real_dist}~=0
%         r_lac_ci(size(r_lac_ci, 1)+1, 1) = cellstr(rogers_cell{ii, idx_lac_ci});
%     end
% end

idx_cell = 2;
idx_lac = 3;
idx_ci = 4;

rural_orangeville_cell = cell(0);
rural_orangeville_ta_radius = [];
rural_orangeville_oh_radius = [];
rural_orangeville_toh_radius = [];
rural_orangeville_isd_radius = [];
rural_orangeville_voro_radius = [];
for ii=1:size(rural_orangeville_tags, 1)
    
    cell_id = rural_orangeville_tags{ii, idx_cell};
    
    cell_split = regexp(cell_id, '-', 'split');
    
    lac_ci = [cell_split{idx_lac}, '-', cell_split{idx_ci}];
    
    loc = find(ismember(t_lac_ci, lac_ci)>0);
    
    rural_orangeville_cell(ii, :) = rogers_cell(loc, :);
    
    loc_radius = find(ismember(r_lac_ci, lac_ci)>0);
    
    rural_orangeville_ta_radius(ii, :) = rogers_ta_radius(loc_radius, :);
    rural_orangeville_oh_radius(ii, :) = rogers_oh_radius(loc_radius, :);
    rural_orangeville_toh_radius(ii, :) = rogers_toh_radius(loc_radius, :);
    rural_orangeville_isd_radius(ii, :) = rogers_isd_radius(loc_radius, :);
    rural_orangeville_voro_radius(ii, :) = rogers_voro_radius(loc_radius, :);
end