function result = inArcDirection(c, cell_angle, p)
    result = 0;
    
    [start_angle, stop_angle] = changeGeoAngle(cell_angle(1), cell_angle(2));
    
%     start_angle = degtorad(90-cell_angle(1));
%     stop_angle = degtorad(90-cell_angle(2));
%     
%     if start_angle<0
%         start_angle = start_angle + 2.*pi;
%     end
%     
%     if stop_angle<0
%         stop_angle = stop_angle + 2.*pi;
%     end
%     
%     % for the angle is clockwise
%     % but, in plane geo is anticlockwise
%     if start_angle<stop_angle
%         start_angle = start_angle + 2.*pi;
%     elseif start_angle==stop_angle
%         result = 1; % it's circle
%         return ;
%     end
    
    geo_angle = [start_angle, stop_angle];
    
    theta = (geo_angle(2) - geo_angle(1))/2;
    
    if theta==pi % it's circle
        result = 1;
        return ;
    end
    
    
    v_p0.x = 1000*cos(sum(geo_angle)/2);
    v_p0.y = 1000*sin(sum(geo_angle)/2);
    
    
    v_p1.x = p.x - c.x;
    v_p1.y = p.y - c.y;
    
    %v_p0.x * v_p1.x + v_p0.y * v_p1.y
    dotP = v_p0.x*v_p1.x + v_p0.y*v_p1.y;
    
    len = sqrt((v_p0.x.^2 + v_p0.y.^2)*(v_p1.x.^2+v_p1.y.^2));
    
    cosValue = dotP/len;
    
    if cosValue<-1
        cosValue = -1;
    elseif cosValue>1
        cosValue = 1;
    end
    
    insec_angle = acos(cosValue);
   
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