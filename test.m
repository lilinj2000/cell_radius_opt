
clc;
% clear;

% figure;
% 
% plot(relative_radius(:, 1:2));
% legend('TA', 'OH');
% 
% figure;
% plot(relative_radius(:, 3:5));
% legend('ISD', 'voronoi_max', 'voronoi_mean');

% pie_legend = {'500', '1000', '3000', '10000', '20000', 'other'};
% M = [250, 750, 1250, 4750, 15250, 24750];
% pie_legend = {'500', '3000', '10000', 'other'};
% M = [250, 750, 5250, 14750];
% N = hist(abs(relative_radius(:, 1)), M);
% 
% pie(N);
% z = N==0;
% 
% pie_index = 1;
% for jj = 1: length(z)
%     if ~z(jj)
%         legend_str{pie_index} = pie_legend{jj};
%         pie_index = pie_index + 1;
%     end
% end
% legend(legend_str, 'Location', 'SouthEast');
% title(methods_name{ii});

% hold on;

% angle = linspace(0, 2.*pi, 1000);
% a = 100;
% b = 50;
% 
% 
% % x = 60 + a.*cos(angle);
% % y = 70 + b.*sin(angle);
% % 
% % plot(x, y, '-');
% 
% plot(60, 70, '+');
% 
% [x, y] = calculateEllipse(60, 70, 100, 50, 60);
% plot(x, y, '-k');
% 
% 
% radius = [ 100, 80, 130];
% line_style = {'k-', 'g--', 'b-.'};
% 
% for ii=1:length(radius)
%     
%     dist = radius(ii);
%     s = line_style{ii};
%     
%     [x, y] = drawCircle(60, 70, dist);
%     plot(x, y, s);
% 
%     deg = rand(1);
%     x1 = dist.*cos(deg.*2.*pi);
%     y1 = dist.*sin(deg.*2.*pi);
%     quiver(60, 70, x1, y1, 0, s);
% end
% 
% min = -100;
% max = 210;
% axis([min, max, min, max]);

% e = [5.0000   43.6487  -79.3687  464.0000  149.0000   31.0000   80.0000];
% 
% figure;
% hold on;
% % rad = deg2rad(e(60));
% %x = linspace(e(2)-e(4), e(2)+e(4), 1000);
% x = linspace(43-464, 43+464, 1000);
% % ezplot('((x-e(1)).^2./(e(4).^2) + (y-e(2)).^2./(e(5).^2))-1');
% ezplot('((x-43).^2./464.^2+(y+79).^2./149.^2)-1', [43-464, 43+464]);
% 
% rad = deg2rad(60);
% % x1 = x.*cos(rad) + y.*sin(rad);
% % y1 = -1.*x.sin(rad) + y.cos(rad);
% 
% % cos_value = cos(rad); %0.5
% % sin_value = sin(rad); %0.866
% 
% ezplot('(((x.*0.5+y.*0.866)-(43.*0.5-79.*0.866)).^2./464.^2 + ((-1.*x.*0.866+y.*0.5)+(-43.*0.866-79.*0.5)).^2./149.^2)-1', [43.*0.5-79.*0.866-464, 43.*0.5-79.*0.866+464]);

% b= (((x.*0.5+y.*0.866)-(43.*0.5-79.*0.866)).^2./464.^2 + ((-1.*x.*0.866+y.*0.5)+(-43.*0.866-79.*0.5)).^2./149.^2)-1;

% b = (((x.*0.5+y.*0.866)-43).^2./464.^2 + ((-1.*x.*0.866+y.*0.5)+79).^2./149.^2)-1

% x = 1:10;
% ezplot('x-y.^2');

% x = gallery('uniformdata',[1 10],0);
% y = gallery('uniformdata',[1 10],1);
% voronoi(x,y)

%pp = rand(6,2).*100;

% axis_max = max( [pp(:,1); pp(:, 2)]);
% axis_min = min( [pp(:,1); pp(:, 2)]);
% 
% % distance = [];
% 
% for ii=1:size(pp, 1)
%     for jj=1:size(pp, 1)
%         distance(ii, jj) = sqrt((pp(ii,1)-pp(jj,1)).^2 + (pp(ii,2)-pp(jj,2)).^2);
%     end
% end
% 
% for ii=1:size(pp,1)
%     isd(ii) = sum(distance(ii, :))./(size(pp,1)-1);
% end
% 
% 
% index = 3;
% 
% figure;
% hold on;
% plot(pp(:, 1), pp(:,2), 'o');
% 
% % draw line
% for ii=1:size(pp, 1)
%     if ii~=index
%         k = (pp(index,2)-pp(ii,2))./(pp(index,1)-pp(ii,1));
%         
%         x = linspace(pp(index,1), pp(ii,1), 1000);
%         y = k.*(x-pp(index,1)) + pp(index,2);
%         
%         plot(x, y, '-k');
%     end
% end
% 
% k = [0.7, 1, 1.2];
% s = {'b-.', 'k-', 'g--'};
% 
% for ii=1:length(k)
% 
% % draw cicle
% angle = 0:2.*pi/100:2.*pi;
% dist = isd(index).*k(ii);
% x = pp(index,1) + dist.*cos(angle);
% y = pp(index,2) + dist.*sin(angle);
% plot(x, y, s{ii});
% 
% % ray
% deg = rand(1);
% x1 = dist.*cos(deg.*2.*pi);
% y1 = dist.*sin(deg.*2.*pi);
% quiver(pp(index,1), pp(index,2), x1, y1, 0, s{ii});
% 
% end
% 
% axis_max = axis_max + isd(index)+10-15;
% axis_min = axis_min - isd(index)-10+15;
% 
% axis([axis_min, axis_max, axis_min, axis_max]);


% angle = linspace(0, 2.*pi, 1000);
% 
% center.x = 3;
% center.y = 3;
% 
% S = cell(3,1);
% S{1} = 'r-';
% S{2} = 'b-';
% S{3} = 'c-';
% 
% deg = randn(3,1);
% 
% figure;
% hold on;
% 
% plot(center.x, center.y, '+');
% 
% index = 0;
% for radius = 5:3:12
%     x = center.x + radius.*cos(angle);
%     y = center.y + radius.*sin(angle);
%     
%     index = index+1;
%     plot(x, y, S{index});
%     
%     % ray
%     x1 = radius.*cos(deg(index).*2.*pi);
%     y1 = radius.*sin(deg(index).*2.*pi);
%     quiver(center.x, center.y, x1, y1, 0, S{index});
%     
% end
% 
% % axis off;
% hold off;
