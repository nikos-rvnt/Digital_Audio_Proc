%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                   DPMS SESE  -  Psifiakh Texnologia Hxou                %
%                                                                         %
%           Ylopoihsh Sympiesth Dynamikhs Perioxhs (Compressor)           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all ;

cd 'audio samples';
d = dir ;
d(1:2) = [] ;

% str = {d.name} ;
% [s,v] = listdlg('PromptString','Epilekste ena hxhtiko deigma:',...
%                 'SelectionMode','single',...
%                 'ListSize',[240 380],'ListString', str) ; 
 
epilogh = d(3).name ;  

% load stereo sound
[hxhtiko,Fs] = audioread(epilogh) ;

% stereo2mono
if( size( hxhtiko, 2)>1)
    hxhtiko = (hxhtiko(:,1) + hxhtiko(:,2))./2 ;
end

cd ..

% normalise input
hxhtiko = normalisedVEC(hxhtiko) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------- Feed-Forward Side Chain ---------------------------

% abs-value 
hxhtikoLD = max( abs(hxhtiko), 1e-6) ;

figure(1)
plot(hxhtiko)
xlabel('samples(n)') ;
ylabel('Amplitude') ;
title('Compressor Normalised Input [-1,1]') ;

figure(2)
plot(hxhtikoLD)
xlabel('samples(n)') ;
ylabel('Amplitude') ;
title('Input Absolute Value') ;

% linear to dB conversion
hxhtiko_dB = 20*log10(hxhtikoLD) ;

figure(3)
plot(hxhtiko_dB)
xlabel('samples(n)') ;
ylabel('Amplitude(dB)') ;
title('After the convertion to dB scale') ;


% choice -> apanthsh tou xrhsth meta apo orismo idiothtwn compressora
choice = '' ;
while( ~strcmp(choice,'Yes') )

    prompt = {'Threshold(dB):','Ratio:','Attack Time(sec.):',...
        'Release Time(sec.):','Knee Width(dB):','Level Detection(rms/peak):'} ;
    dlg_title = 'Compressor' ;
    num_lines = 1 ;
    defaultans = {'-25','10','0.0015','0.0180','10','peak'} ;
    options.Resize = 'on' ;
    options.WindowStyle = 'normal' ;
    size_wind = [1 50; 1 50; 1 50; 1 50] ;
    xrhsths = inputdlg( prompt, dlg_title, num_lines, defaultans, options) ;

    Threshold = str2num(xrhsths{1,1}) ;
    Ratio = str2num(xrhsths{2,1}) ;
    T_att = str2num(xrhsths{3,1}) ;
    T_rel = str2num(xrhsths{4,1}) ;
    KneeW = str2num(xrhsths{5,1}) ;
    LVLD =  (xrhsths{6,1}) ;

% -------------------------- Gain Computer --------------------------------

% Threshold: 
%Threshold = -45 ;

% Ratio: 
%Ratio = 8 ;

% Knee Width:
%KneeW = 10 ;

    inputEX = -60:5:0 ;
    inputEX = inputEX' ; 
    inputG = hxhtiko_dB ; 
    
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
        % case3: 2*(xG - T) > W
        elseif( 2*(inputEX(i,1)-Threshold) > KneeW )
            outputEX(i,1) = Threshold + (inputEX(i,1) - Threshold)/Ratio ;    
        end
    
    end
end

   figure(4)
   plot( inputEX, outputEX)
   hold on ;
   plot( inputEX, inputEX)
   plot( Threshold*ones( size(inputEX,2), size(inputEX,1)),inputEX)
   hold off ;
   xlabel('Input (dB)') ;
   ylabel('Output (dB)') ;
   legend(['Compression with Ratio=',num2str(Ratio),':1'],'No Compression(1:1)','Thhreshold') ;
   title('Compressor Static Characteristic Curve') ;

     options.Interpreter = 'tex';
     % Include the desired Default answer
     options.Default = 'Yes';
     % Use the TeX interpreter in the question
     qstring = 'Do you want to apply the Compressor with these settings ?';
     choice = questdlg( qstring,'Apply Compressor',...
          'Yes','Change Settings', options) ;

end 

inputG = hxhtiko_dB ; 

if(isequal(KneeW,0))
    % hard knee characteristic
    outputG(inputG<=Threshold) = inputG(inputG<=Threshold) ;
    outputG(inputG>Threshold) = Threshold + (inputG(inputG>Threshold) - Threshold)/Ratio ;

else
    % soft knee characteristic
    for(i=1:length(inputG'))
    
        % case1 : 2*(xG - T) < - W
        if( 2*(inputG(i,1)-Threshold) < -KneeW )
            outputG(i,1) = inputG(i,1) ;
        % case2: 2*abs(xG-T) <= W    
        elseif( 2*abs(inputG(i,1)-Threshold) <= KneeW )
            outputG(i,1) = inputG(i,1)+(((1/Ratio)-1)*((inputG(i,1)-Threshold+KneeW/2)^2))/(2*KneeW) ;
        % case3: 2*(xG - T) > W
        elseif( 2*(inputG(i,1)-Threshold) > KneeW )
            outputG(i,1) = Threshold + (inputG(i,1) - Threshold)/Ratio ;    
        end
     
     coga(i,1) = -outputG(i,1) + inputG(i,1) ;
    
    end   
end

figure(5)
plot(-coga)
xlabel('samples(n)') ;
ylabel('Amplitude(dB)') ;
title('Computed Gain') ;


% --------------------------- gain smoothing ------------------------------

% attack time coefficient
Co_att = exp((-log(9))/(Fs*T_att)) ;
% release time coefficient
Co_rel = exp((-log(9))/(Fs*T_rel)) ;


% level detection - gain smoothing
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
    smga( 1, 1) = sqrt( coga( 1, 1).^2 ) ;
        
    for i=2:length(coga)
        % square of the level is an estimate of square of the input
        smga(i,1) = sqrt( Co_att*smga( i-1, 1)^2 + (1-Co_att)*coga( i, 1)^2 ) ;        
    end
        
end


figure(6)
subplot(2,1,1)
plot(-coga)
xlabel('samples(n)') ;
ylabel('Amplitude(dB)') ;
title('Before the Level Detection/Gain Smoothing') ;
subplot(2,1,2)
plot(-smga)
xlabel('samples(n)') ;
ylabel('Amplitude(dB)') ;
title('After the Level Detection/Gain Smoothing') ;

%% -------------------- Make-Up Gain Addition ------------------------------

% auto make up gain computation - negative of the computed gain for 0dB

if(KneeW/2<Threshold)
    tempMUG = 0 ;
    
elseif(2*abs(Threshold)<=KneeW)
    tempMUG = (1/Ratio - 1)*(Threshold-(W/2))^2 ;

elseif(Threshold<(-KneeW/2))
    tempMUG = Threshold + Threshold/Ratio ;
    
end

% make-up gain
magain = -smga  -tempMUG ; % muGainS ;

figure(7)
subplot(2,1,1)
plot(-smga)
xlabel('samples(n)') ;
ylabel('Amplitude(dB)') ;
title('Before Auto Make-Up Gain') ;
subplot(2,1,2)
plot(magain)
xlabel('samples(n)') ;
ylabel('Amplitude(dB)') ;
title('After to Auto Make-Up Gain') ;

% --------------------------- dB 2 Linear ---------------------------------
lingain = 10.^(magain/20) ;

ThresLin = 10.^(Threshold/20) ;

% ------------------------ Compressed Output ------------------------------
compOUT = hxhtiko.*lingain ;

% normalise output
compOUTn = normalisedVEC(compOUT) ;


figure(8)
subplot(1,2,1)
plot(hxhtiko)
hold on ;
plot( ThresLin*ones( size(compOUT) ) )
hold off ;
title('Compressor Input') ;
xlabel('Samples(n)') ;
ylabel('Amplitude') ;
legend('Input','Threshold')
subplot(1,2,2)
plot(compOUTn)
hold on ;
plot( ThresLin*ones( size(compOUT) ) )
hold off ;
title('Compressor Output') ;
xlabel('Samples(n)') ;
ylabel('Amplitude') ;
legend('Output','Threshold') ;
