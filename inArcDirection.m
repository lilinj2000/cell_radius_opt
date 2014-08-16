function result = inArcDirection(c, cell_angle, p)
    result = 0;
    
    start_angle = degtorad(90-cell_angle(1));
    stop_angle = degtorad(90-cell_angle(2));
    
    if start_angle<0
        start_angle = start_angle + 2.*pi;
    end
    
    if stop_angle<0
        stop_angle = stop_angle + 2.*pi;
    end
    
    % for the angle is clockwise
    % but, in plane geo is anticlockwise
    if start_angle<stop_angle
        start_angle = start_angle + 2.*pi;
    end
    
    geo_angle = [stop_angle, start_angle];
    
    theta = (geo_angle(2) - geo_angle(1))/2;
    
    if theta >= pi % it's circle
        result = 1;
        return ;
    end
    
    p0.x = c.x + 1000*cos(sum(geo_angle)/2);
    p0.y = c.y + 1000*sin(sum(geo_angle)/2);
    
    
    %p.x * p0.x + p.y * p0.y
    dotP = p.x*p0.x + p.y*p0.y;
    len = sqrt((p.x.^2 + p.y.^2)*(p0.x.^2+p0.y.^2));
    
    insec_angle = acos(dotP/len);
   
    if insec_angle<=theta
        result = 1;
    end
    
%     figure;
%     hold on;
%     drawArc(c, 1000, geo_angle(1), geo_angle(2));
%     
%     drawLine(c, p0);
%     drawLine(c, p);
%     plot(p.x, p.y, '*');
%     axis equal;    
    
end