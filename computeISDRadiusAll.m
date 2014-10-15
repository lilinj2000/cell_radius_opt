%% compute isd radius all

function isd_radius_all = computeISDRadiusAll(cell_enu)

idx_lat = 1;
idx_long = 2;


log = log4m.getLogger('isdall.log');
log.setCommandWindowLevel(log.OFF);

isd_radius_all = zeros(size(cell_enu, 1));

for ii = 1 : size(cell_enu, 1)
%     cell_enu_angle_dist = [];
    
    log.info('computeISDRadiusAll()', ['compute no ', num2str(ii)]);

    for jj = 1 : size(cell_enu, 1)
        if ii~=jj
            if isd_radius_all(jj,ii)~=0
                isd_radius_all(ii, jj)= isd_radius_all(jj, ii);
            else
                isd_radius_all(ii, jj) = sqrt((cell_enu(ii, idx_lat)-cell_enu(jj, idx_lat)).^2 ...
                    + (cell_enu(ii, idx_long)-cell_enu(jj, idx_long)).^2);
            end
        end
    end
end

end