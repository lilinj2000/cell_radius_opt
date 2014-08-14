function similar_value = similarAlgorithm(data1, data2 )
% input two sets of data, calculate the similar value.
% 
% similar algorithm:
%        i from 1 to  n
%  dx = sqrt(sum((data1(i) - data2(i)).^2))
% similar_value = 1./(1+dx);

similar_value = 0;
if length(data1) ~= length(data2)
    return ;
end

%
% sum = 0;
% for ii = 1: length(data1)
%     sum = sum + (data1(ii)-data2(ii)).^2;
% end
% 
% dx = sqrt(sum);
% 
% similar_value = 1./(1+dx);


% Tanimoto Coefficient
% a = sum(xy);
% b = sum(x.^2);
% c = sum(y.^2);
% similar = a./(b+c-a);

a = sum(data1.*data2);
b = sum(data1.^2);
c = sum(data2.^2);

similar_value = a./(b+c-a);

end