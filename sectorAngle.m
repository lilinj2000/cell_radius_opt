function [start_angle, stop_angle] = sectorAngle(antenna_type, cell_dir, sector_angle)

% angle process
if strcmp(antenna_type, 'OMNI')
    start_angle = 0;
    stop_angle = 360;
else
    
    if sector_angle==0
        sector_angle = 120;
    end

    if cell_dir<sector_angle./2
        start_angle = 360 - (sector_angle./2-cell_dir);
    else
        start_angle = cell_dir - sector_angle./2;
    end

    stop_angle = cell_dir + sector_angle./2;
end
    
end