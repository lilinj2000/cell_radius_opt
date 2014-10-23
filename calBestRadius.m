function [radius, similar_best, k] = calBestRadius(geo_radius, r_radius)
% calculate the best radius according the geo radius, and the real radius
% the constraint is the similary best

steps = [0.8:0.1:5];

radius = [];
similar_best = 0;
k = [];

for step=steps
    tmp_radius = geo_radius*step;
    
    similar_value = similarAlgorithm(tmp_radius, r_radius);
    
    if similar_value>similar_best
        similar_best = similar_value;
        radius = tmp_radius;
        k = step;
    end
end

end