function data = freqToNumber(cell_freq)
    switch cell_freq
        case 'GSM800'
            data = 10;
        case 'GSM900'
            data = 20;
        case 'GSM1800'
            data = 30;
        case 'GSM1900'
            data = 40;
        otherwise
            data = 10;
    end
end