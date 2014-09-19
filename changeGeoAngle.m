function [dst_start_angle, dst_stop_angle] = changeGeoAngle(src_start_angle, src_stop_angle)

deg_start_angle = degtorad(90-src_start_angle);
    
if deg_start_angle<0
    deg_start_angle = deg_start_angle + 2.*pi;
end

deg_stop_angle = degtorad(90-src_stop_angle);
    
if deg_stop_angle<0
    deg_stop_angle = deg_stop_angle + 2.*pi;
end

if deg_start_angle<=deg_stop_angle
    deg_start_angle = deg_start_angle + 2.*pi;
end

% for the angle is clockwise
% but, in plane geo is anticlockwise
[dst_start_angle, dst_stop_angle] = deal(deg_stop_angle, deg_start_angle);


end