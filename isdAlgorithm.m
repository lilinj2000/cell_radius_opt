function radius_outer = isdAlgorithm(cell_enu, neighbour_count, isd_clearance)

radius_outer = zeros(length(cell_enu), 1); 

distance = zeros(length(cell_enu));

for ii = 1 : length(cell_enu)
    for jj = 1 : length(cell_enu)
        if ii~=jj
            distance(ii, jj) = sqrt((cell_enu(ii, 1)-cell_enu(jj, 1)).^2 ...
                + (cell_enu(ii, 2)-cell_enu(jj, 2)).^2);
        end
    end
    
    [~, index] = sort(distance(ii, :));
    
    % get rif off the distance which is less than isd_clearance
    index_small = length(index)+1;
    for jj = 1 : length(index)
        if distance(ii, index(jj))>isd_clearance
            index_small = jj;
            break;
        end
    end
    
    count = length(index)-index_small+1;
    if count>neighbour_count
        count = neighbour_count;
    end
    
    radius_outer(ii) = sum(distance(ii, index(index_small:index_small+count-1))) ...
        ./count;
        
end

end