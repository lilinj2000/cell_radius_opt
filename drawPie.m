function result = drawPie(data, M, pie_legend, pie_title)

    if length(M)~=length(pie_legend)
        error('invalid input');
    end

    
    N = hist(data, M);

    figure;
    pie(N);
    
    
    z = N==0;

    pie_index = 1;
    for jj = 1: length(z)
        if ~z(jj)
            legend_str{pie_index} = pie_legend{jj};
            pie_index = pie_index + 1;
        end
    end
    
    legend(legend_str, 'Location','EastOutside');
    
    title(pie_title);
    
    result = round(N./sum(N).*100) ;
    
end