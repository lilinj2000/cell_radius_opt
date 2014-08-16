% regression analysis

X = [];
% for ii=1:length(cell_data)
% %     X(ii, 1) = log10(freqToNumber(char(cell_data{ii, CELL_FREQ})));
% %     X(ii, 2) = (log10(freqToNumber(char(cell_data{ii, CELL_FREQ})))).^2;
% %     X(ii, 3) = exp(classToNumber(char(cell_data{ii, CELL_CLASS})));
%     X(ii, 1) = power(10, (str2double(cell_data{ii, CELL_POWER}) ...
%         + str2double(cell_data{ii, CELL_ACCMIN})+21-116)/35);
% end

% X(:, 1) = log(radius(:, 4));
% X(:, 2) = (log(radius(:, 4))).^2;
% X(:, 3) = radius(:, 4);
X = radius(:, 3:6);

% X(:, 2) = exp(radius(:,5));
% X(:, 2) = log(radius(:,5));


total_data = length(real_radius);
stat_section = [50, 100, 200, 300, 400, 500, 1000, 2000, 3000, 4000, 5000];
ration_data = [];


% Multiple linear regression
index_method = 1;
[b, bint, r, rint, stats] = regress(real_radius', ...
    [ones(length(real_radius), 1) X], 0.2);
dst_data = b(1) + b(2).*X(:, 1) + b(3).*X(:, 2) + b(4).*X(:, 3)+ b(5).*X(:, 4);
relative_data = dst_data - real_radius';

for ii=1:length(stat_section)
    ration_data(index_method, ii) = length(find(abs(relative_data)<=stat_section(ii)))/total_data*100;
end


 
% %Polynomial curve fitting
% index_method = 2;
% [p, S] = polyfit(observe_data, real_data, 1);
% dst_data = polyval(p, observe_data);
% relative_data = dst_data - real_data;
% for ii=1:length(stat_section)
%     ration_data(index_method, ii) = length(find(abs(relative_data)<=stat_section(ii)))/total_data*100;
% end
% 
% %Polynomial curve fitting
% index_method = 3;
% [p, S] = polyfit(observe_data, real_data, 2);
% dst_data = polyval(p, observe_data);
% relative_data = dst_data - real_data;
% for ii=1:length(stat_section)
%     ration_data(index_method, ii) = length(find(abs(relative_data)<=stat_section(ii)))/total_data*100;
% end
% 
% %Polynomial curve fitting
% index_method = 4;
% [p, S] = polyfit(observe_data, real_data, 3);
% dst_data = polyval(p, observe_data);
% relative_data = dst_data - real_data;
% for ii=1:length(stat_section)
%     ration_data(index_method, ii) = length(find(abs(relative_data)<=stat_section(ii)))/total_data*100;
% end