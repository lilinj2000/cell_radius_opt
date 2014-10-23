% cal the best radius

idx_radius = 3;

% beijing
% load beijing_radius;
% isd_radius = beijing_isd_radius;
% voro_radius = beijing_voro_radius;
% real_radius = r_real_radius(:, idx_radius);
% 

% turkey
% load turkey_radius;
% isd_radius = turkey_isd_radius_95;
% voro_radius = turkey_voro_radius_95;
% real_radius = r_real_radius(:, idx_radius);
% 

% rogers
% load rogers_radius
isd_radius = rogers_isd_radius_95;
voro_radius = rogers_voro_radius_95;
real_radius = r_real_radius(:, idx_radius);



isd_best_radius = [];
isd_similar_best = [];
isd_k = [];
for ii=1:size(isd_radius, 2)
    [isd_best_radius(:, ii), isd_similar_best(:, ii), isd_k(:, ii)] = ...
        calBestRadius(isd_radius(:, ii), real_radius);
end

voro_best_radius = [];
voro_similar_best = [];
voro_k = [];
for ii=1:size(voro_radius, 2)
    [voro_best_radius(:, ii), voro_similar_best(:, ii), voro_k(:, ii)] = ...
        calBestRadius(voro_radius(:, ii), real_radius);
end
