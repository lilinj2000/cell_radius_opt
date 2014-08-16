function data = classToNumber(cell_class)
    switch cell_class
       case 'MACRO'
            data = 30;
       case 'MICRO'
            data = 20;
        otherwise
            data = 10;
    end
end