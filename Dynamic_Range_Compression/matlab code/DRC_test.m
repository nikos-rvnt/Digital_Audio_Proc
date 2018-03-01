function [ compOUTS, smga, coga] = DRC_test(input, Threshold, KneeW, Ratio, Tatt, Trel,LVLD,Fs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                   DPMS SESE  -  Psifiakh Texnologia Hxou                %
%                                                                         %
%           Ylopoihsh Sympiesth Dynamikhs Perioxhs (Compressor)           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --------------------- Feed-Forward Side Chain ---------------------------

Threshold = Threshold ; % str2num(xrhsths{1,1}) ;
Ratio = Ratio ; %4 ; % str2num(xrhsths{2,1}) ;
T_att = Tatt ; %0.8 ; % str2num(xrhsths{3,1}) ;
T_rel = Trel ; % 0.34 ; % str2num(xrhsths{4,1}) ;
KneeW = KneeW ; % 8 ; % str2num(xrhsths{5,1}) ;

% level detection 
hxhtikoLD = max( abs(input), 1e-6) ; 

% linear to dB conversion
hxhtiko_dB = 20*log10(hxhtikoLD) ; 


% -------------------------- Gain Computer --------------------------------

% Threshold 
%Threshold = -45 ;

% Ratio 
%Ratio = 8 ;

% Knee Width 
%KneeW = 10 ;

inputEX = -60:5:0 ;
inputEX = inputEX' ; 
    
if(isequal(KneeW,0))
    outputEX(inputEX<=Threshold) = inputEX(inputEX<=Threshold) ;
    outputEX(inputEX>Threshold) = Threshold + (inputEX(inputEX>Threshold) - Threshold)/Ratio ;

else

    for(i=1:length(inputEX))

        % case1 : 2*(xG - T) < - W
        if( 2*(inputEX(i,1)-Threshold) < -KneeW )
            outputEX(i,1) = inputEX(i,1) ;
        % case2: 2*abs(xG-T) <= W    
        elseif( 2*abs(inputEX(i,1)-Threshold) <= KneeW )
            outputEX(i,1) = inputEX(i,1)+(((1/Ratio)-1)*((inputEX(i,1)-Threshold+KneeW/2)^2))/(2*KneeW) ;
        % case3: 2*(xG - T) > -W
        elseif( 2*(inputEX(i,1)-Threshold) > KneeW )
            outputEX(i,1) = Threshold + (inputEX(i,1) - Threshold)/Ratio ;    
        end

    end

end

figure(1)
plot( inputEX, outputEX)
hold on ;
plot( inputEX, inputEX)
plot( Threshold*ones( size(inputEX,2), size(inputEX,1)),inputEX)
hold off ;
xlabel('Input (dB)') ;
ylabel('Output (dB)') ;
legend(['Compression with Ratio=',num2str(Ratio),':1'],'No Compression(1:1)','Thhreshold') ;
title('Compressor Static Characteristic Curve') ;


inputG = hxhtiko_dB' ; 

% hard knee characteristic - hard half-wave rectifier
if(isequal(KneeW,0))
    
    for(i=1:length(inputG))
    
        if(inputG(i,1)<=Threshold)
        
            outputG(i,1) = inputG(i,1) ;
        else
            
            outputG(i,1) = Threshold + (inputG(i,1) - Threshold)/Ratio ;      
        end
             
        % computer gain on soft knee characteristic
        coga(i,1) = -outputG(i,1) + inputG(i,1) ;
    end
    
else
            
     for(i=1:length(inputG))
            
        % soft knee characteristic - soft half-wave rectifier
    
        % case1 : 2*(xG - T) < - W
        if( 2*(inputG(i,1)-Threshold) < -KneeW )
            outputG(i,1) = inputG(i,1) ;
        % case2: 2*abs(xG-T) <= W    
        elseif(( 2*abs(inputG(i,1)-Threshold) <= KneeW ))
            outputG(i,1) = inputG(i,1)+((((1/Ratio)-1)*((inputG(i,1)-Threshold+(KneeW/2))^2))/(2*KneeW)) ;
        % case3: 2*(xG - T) > W
        elseif( 2*(inputG(i,1)-Threshold) > KneeW )
            outputG(i,1) = (Threshold + (inputG(i,1) - Threshold)/Ratio) ;    
        end

        % computer gain on soft knee characteristic
        coga(i,1) = -outputG(i,1) + inputG(i,1) ;
    end    
  
end


%% --------------------------- gain smoothing ------------------------------

% attack time coefficient
Co_att = exp((-log(9))/(Fs*T_att)) 
% release time coefficient
Co_rel = exp((-log(9))/(Fs*T_rel)) 

smga = zeros( size(coga) ) ;
% level detection 
if( strcmp(LVLD,'peak') )

    if(coga( 1, 1)>0)
        smga( 1, 1) = (1 - Co_att)*coga( 1, 1) ;

    else
        smga( 1, 1) = (1 - Co_rel)*coga( 1, 1) ;
    
    end
    
    % smooth branching peak detector
    for (i=2:length(coga))

        % one pole-low pass filter
        if (coga(i,1) > smga( i-1, 1))
            smga(i,1) = Co_att*smga( i-1, 1) + (1 - Co_att)*(coga(i)) ;
        else
            smga(i,1) = Co_rel*smga( i-1, 1) + (1 - Co_rel)*(coga(i)) ;
        end

    end
    
elseif( strcmp(LVLD,'rms') )
    % rms based level detection 
    smga( 1, 1) = sqrt( coga( 1, 1)^2 ) ;
        
    for i=2:length(coga)
        smga(i,1) = sqrt( Co_att*smga(i-1,1)^2 + (1-Co_att)*coga(i,1)^2 ) ; 
        
    end
        
end


% -------------------- Make-Up Gain Addition ------------------------------

% auto make up gain computation - negative of the computed gain for 0dB

if(KneeW/2<Threshold)
    tempMUG = 0 ;
    
elseif(2*abs(Threshold)<=KneeW)
    tempMUG = (1/Ratio - 1)*(Threshold-(KneeW/2))^2 ;

elseif(Threshold<(-KneeW/2))
    tempMUG = Threshold + Threshold/Ratio ;
    
end

magain = -smga - tempMUG ; % muGainS ;

%magain = -smga + 10 ;

% --------------------------- dB 2 Linear --------------------------------- 

lingain = 10.^(magain/20) ;

ThresLin = 10.^(Threshold/20) ;

% ------------------------ Compressed Output ------------------------------
compOUTS = input.*(lingain') ;

% normalise output
compOUTSn = normalisedVEC(compOUTS) ;
inputN = normalisedVEC(input) ;


end
