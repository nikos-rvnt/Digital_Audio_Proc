function dataOUT = normalisedVEC( dataIN )
    % synarthsh kanonikopoihshs timwn dianysmatos sto diasthma [-1,1]
    % diathreitai to arxiko proshmo twn timwn tou dianysmatos

    if(abs( min(dataIN) ) > max(dataIN) )

        megisth_timh = abs( min(dataIN) ) ;
        elaxisth_timh = min(dataIN) ;
    else
        megisth_timh = max(dataIN) ;
        elaxisth_timh = -max(dataIN) ;
    
    end

    % xout = 2*(xin - min)/(max - min) - 1 -> (b-a)*(xin - min)/(max -
    %           min) + a
    
    dataOUT = 2.*(dataIN - elaxisth_timh)./(megisth_timh - elaxisth_timh) - 1 ;

end
