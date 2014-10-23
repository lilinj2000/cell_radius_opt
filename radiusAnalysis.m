function radiusAnalysis(radius, cell_lat_long, subhead)

figure; % plot the cell geo location
hold on;

radius_section = [0, 500, 1000, 2000, 4000, 6000, 8000, 10000];
colors = {'b^', 'g^', 'r^', 'b*', 'g*', 'r*', 'b+', 'g+'};
legend_string = {'500', '1000', '2000', '4000', '6000', '8000', ...
    '10000', '> 10000'};
legend_h = zeros(1,length(legend_string));

% plot(cell_lat_long(:, 1), cell_lat_long(:, 2), '^');
for ii = 1:length(cell_lat_long)
    
    for jj=length(radius_section) : -1 : 1
        if radius(ii)> radius_section(jj)
            h = plot(cell_lat_long(ii, 1), cell_lat_long(ii, 2), colors{jj});
            
            if legend_h(jj)==0
                legend_h(jj) = h;
            end
            break;
        end
    end
end

real_legend_h = [];
real_legend_string = cell(0);
for ii = 1:length(legend_h)
    if legend_h(ii)~=0
        index = length(real_legend_h) + 1;
        real_legend_h(index) = legend_h(ii);
        real_legend_string{index} = legend_string{ii};
    end
end
legend(real_legend_h, real_legend_string, 'Location','EastOutside');
title(vertcat({'geo-location'}, subhead));

figure; % radius scatter dia
% plot the radius
plot(radius, '*');
grid on;
title(vertcat({'radius scatter'}, subhead));

figure; % radius bar dia
% 500, 1000, 2000, 4000, 6000, 8000, 10000, 12000
x_lable = {'500', '1000', '2000', '4000', '6000', '8000', '10000', '>10000'};
nbins = [250, 750, 1250, 2750, 5250, 6750, 9250, 10750];
N = hist(radius, nbins);
bar([1 : length(nbins)], N);
for ii = 1 : length(N)
    if N(ii)~=0
        text(ii, N(ii)+3, num2str(N(ii)));
    end
end
set(gca, 'xtick', [1 : length(nbins)]);
set(gca, 'xticklabel', x_lable);
title(vertcat({'radius histogram'}, subhead));
ylim([0, max(N)*1.25]);



% pie 
pie_legend = x_lable ;
pie_title = vertcat({'radius pie'}, subhead);
drawPie(radius, nbins, pie_legend, pie_title);

figure; % cdfplot
cdfplot(radius);
title(vertcat({'radius cdf'}, subhead));


end