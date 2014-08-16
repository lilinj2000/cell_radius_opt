clc;
clear;

close all;

global LAC_CI;
global LAT_LONG;
global CELL_ANGLE;
global CELL_FREQ;
global CELL_CLASS;
global CELL_POWER;
global CUSTOM_RADIUS;
global CELL_TALIM;
global CELL_ACCMIN;
global REAL_RADIUS;


rogers = 1;
turkey = 2;
india = 3;

rogers_area_rural_orangeville = 1;
rogers_area_suburban_brampton = 2;
rogers_area_suburban_scarborough = 3;
rogers_area_suburban = 4;
rogers_area_urban_dt = 5;
rogers_area_all = 6;

subhead = {'(rogers urban dt - cell geo)'};

% init the cell
cell_data = initialCellData(rogers, rogers_area_urban_dt);


for ii=1:length(cell_data)
    cell_enu(ii, 1) = cell_data{ii, LAT_LONG}(1);
    cell_enu(ii, 2) = cell_data{ii, LAT_LONG}(2);
    
end

v = VoronoiFortuneAlgo(cell_enu, 0.5);

v.do();

voronoi_mean_distance = [];

MAX_WAY = 1;
MEAN_WAY = 2;
MIN_WAY = 3;

for ii = 1: length(cell_enu)
    p.x = cell_enu(ii, 1);
    p.y = cell_enu(ii, 2);

    voronoi_mean_distance(ii, 1) = v.calculateRadius(p, MEAN_WAY, cell_data{ii, CELL_ANGLE});
end

radiusAnalysis(voronoi_mean_distance, cell_enu, subhead);
