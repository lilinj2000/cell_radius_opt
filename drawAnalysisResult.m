
function drawAnalysisResult(radius, real_radius, methods_name, subhead)

pie_legend = {'10%', '20%', '30%', '50%', '>50%'};
M = [0.05, 0.15, 0.25, 0.35, 0.65];


pie_legend2 = {'500', '1000', '3000', '6000', '10000', '>10000'};
M2 = [250, 750, 1250, 4750, 7250, 12750];
        
% switch area_type
%     case urban
%         pie_legend2 = {'50', '100', '200', '300', '500', '>500'};
%         M2 = [25, 75, 125, 275, 325, 675];
%         
%     case {suburban, rural, all}
%         pie_legend2 = {'500', '1000', '3000', '6000', '10000', '>10000'};
%         M2 = [250, 750, 1250, 4750, 7250, 12750];
%         
% %         pie_legend2 = {'500', '1000', '2000', '4000', '6000', '8000', '10000', '>10000'};
% %         M2 = [250, 750, 1250, 2750, 5250, 6750, 9250, 10750];
%         
%         % pie_legend = {'500', '1000', '3000', '10000', '20000', 'other'};
%         % M = [250, 750, 1250, 4750, 15250, 24750];
% %         pie_legend2 = {'500', '3000', '10000', 'other'};
% %         M2 = [250, 750, 5250, 14750];
% end

abs_radius_compare = [];
per_radius_compare = [];

for ii = 1 : size(radius, 2)
    pie_title = vertcat(methods_name{ii}, subhead);
    
    relative_radius(:, ii) = radius(:, ii) - real_radius;
    
    abs_radius_compare(:, ii) = drawPie(abs(relative_radius(:, ii)), M2, pie_legend2, pie_title);
    
%     figure;
%     plot(relative_radius(:, ii));
%     legend(methods_name{ii});
%     
%     figure;
%     
%     hist(abs(relative_radius(:, ii)));
%     legend(methods_name{ii});
%     

    relative_radius_percent = abs(relative_radius(:, ii))./real_radius;
    
    
    per_radius_compare(:, ii) = drawPie(relative_radius_percent, M, pie_legend, pie_title);

%     drawPie(radius(:, ii), M, pie_legend, methods_name{ii});
    
%     figure;
%     plot(real_radius, '*');
%     hold on;
%     plot(radius(:, ii), '^');
%     
end


% clear cell;
% radius_by_cell = cell(5,3);
% for ii = 1 : length(cell_chr)
%     row_count = size(radius_by_cell{cell_chr(ii,1), cell_chr(ii,2)}, 1);
%     radius_by_cell{cell_chr(ii,1), cell_chr(ii,2)}(row_count+1, :) = radius(ii, :);
% end
% 
% relative_radius_by_cell = cell(5,3);
% for ii = 1 : size(radius_by_cell, 1)
%     for jj = 1 : size(radius_by_cell, 2)
%         col_count = size(radius_by_cell{ii, jj}, 2);
%         
%         if col_count>0 
%             for kk = 1 : col_count-1
%                 relative_radius_by_cell{ii, jj}(:, kk) = radius_by_cell{ii,jj}(:, kk) ...
%                    - radius_by_cell{ii, jj}(:, col_count); 
%             end
%         end
%     end
% end

% caluculate the similar value
similar = zeros(size(radius, 2), 1);
for ii = 1 : size(radius, 2)
    similar(ii) = similarAlgorithm(real_radius, radius(:, ii));
end

% the real readius
% methods = methods + 1;
% index_tags = methods;
% methods_name{index_tags} = 'real radius';
% 
% radius(:, index_tags) = real_radius;

figure;
plot(similar, 'o--');
% for ii = 1: length(similar)
%     text(ii+0.1, similar(ii), methods_name{ii});
% end
axis([0, length(similar)+1, 0, 1]);
grid on;
set(gca, 'xtick', [1 : length(similar)]);
set(gca, 'xticklabel', methods_name);
title(subhead);



%check result
figure;
bar(abs_radius_compare);
set(gca, 'xtick', [1 : size(abs_radius_compare, 1)]);
set(gca, 'xticklabel', pie_legend2);
title({'radius compare'});
legend(methods_name, 'Location','EastOutside');
grid on;

% per_radius_compare'
figure;
bar(per_radius_compare);
set(gca, 'xtick', [1 : size(per_radius_compare, 1)]);
set(gca, 'xticklabel', pie_legend);
title({'radius compare'});
legend(methods_name, 'Location','EastOutside');
grid on;


abs_radius_compare'
per_radius_compare'
similar

% similar
% b'
% inmodel
% stats

% legend show;


% output the stat info
% figure;
% bar(stat_report);

% output cell info
% figure;
% 
% plot(relative_radius(:, index_ta_radius:index_oh_propagation_model));
% legend(methods_name{index_ta_radius}, methods_name{index_oh_propagation_model});

% figure;
% plot(relative_radius(:, index_isd_radius:index_voronoi_radius_mean));
% legend(methods_name{index_isd_radius}, methods_name{index_voronoi_radius_max}, methods_name{index_voronoi_radius_mean});

% X = radius(:, index_isd_radius:index_voronoi_radius_min);
% stepwise(X, real_radius');

% figure;
% hold on;
% for ii = 1 : length(cell_enu)
% %ii = 2;
%     % draw start line
%     start_angle = degtorad(90-cell_data{ii, CELL_ANGLE}(1));
%     stop_angle = degtorad(90-cell_data{ii, CELL_ANGLE}(2));
%     
%     if start_angle<0
%         start_angle = start_angle + 2.*pi;
%     end
%     
%     if stop_angle<0
%         stop_angle = stop_angle + 2.*pi;
%     end
%     
%     if start_angle<stop_angle
%         start_angle = start_angle + 2.*pi;
%     end
%     
%     x = linspace(cell_enu(ii, 1), cell_enu(ii, 1) + radius_outer(ii)*cos(start_angle), 1000);
%     y = tan(start_angle)*(x-cell_enu(ii, 1)) + cell_enu(ii, 2);
%     plot(x, y);
%     
%     % draw the stop line
%     x = linspace(cell_enu(ii, 1), cell_enu(ii, 1) + radius_outer(ii)*cos(stop_angle), 1000);
%     y = tan(stop_angle)*(x-cell_enu(ii, 1)) + cell_enu(ii, 2);
%     plot(x, y);
%     
%     % draw the center line
%     x = linspace(cell_enu(ii, 1), cell_enu(ii, 1) + radius_outer(ii)*cos((stop_angle+start_angle)/2), 1000);
%     y = tan((stop_angle+start_angle)/2)*(x-cell_enu(ii, 1)) + cell_enu(ii, 2);
%     plot(x, y);
%     
%     
%     % draw arc - outer radius
%     angle = linspace(stop_angle, start_angle, 1000);
%     plot(cell_enu(ii, 1) + radius_outer(ii).*cos(angle), cell_enu(ii, 2) + radius_outer(ii).*sin(angle));
%     
%     
%     % draw arc - inner radius
%     angle = linspace(stop_angle, start_angle, 1000);
%     plot(cell_enu(ii, 1) + radius_inner(ii).*cos(angle), cell_enu(ii, 2) + radius_inner(ii).*sin(angle));
%     
%     break;
% end
% hold off;
% axis equal;

end