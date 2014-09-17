clc;
clear;

close all;

importSysVar;

subhead = {'(china beijing - cell geo)'};

% init the cell
cell_data = initialCellData(china, china_beijing);


for ii=1:length(cell_data)
    cell_enu(ii, 1) = cell_data{ii, LAT_LONG}(1);
    cell_enu(ii, 2) = cell_data{ii, LAT_LONG}(2);
end

v = VoronoiFortuneAlgo(cell_enu, 0.5);

v.do();

voronoi_mean_distance = [];

cell_enu_tag =[];
voronoi_mean_distance = [];
for ii = 1: length(cell_data)
    
    if cell_data{ii, REAL_RADIUS}~=0
        cell_enu_tag = [cell_enu_tag; cell_enu(ii, :)];
        
        p.x = cell_enu(ii, 1);
        p.y = cell_enu(ii, 2);

        [~, mean_distance] = v.calculateRadius(p, cell_data{ii, CELL_ANGLE});
        
        voronoi_mean_distance = [voronoi_mean_distance; mean_distance];
    end
end

radiusAnalysis(voronoi_mean_distance, cell_enu_tag, subhead);
