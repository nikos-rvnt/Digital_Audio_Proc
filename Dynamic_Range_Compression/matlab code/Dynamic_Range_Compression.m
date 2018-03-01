%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                   DPMS SESE  -  Psifiakh Texnologia Hxou                %
%                                                                         %
%           Ylopoihsh Sympiesth Dynamikhs Perioxhs (Compressor)           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


choiceC = 'Again' ;
while(strcmp(choiceC,'Again'))

clear all ;
   
    
cd 'audio samples';
d = dir ;
d(1:2) = [] ;

str = {d.name} ;
[s,v] = listdlg('Name','Compressor with Auto Make-Up Gain',...
                'PromptString',...
                'Choose an audio sample to apply on it the Compressor:',...
                'SelectionMode','single',...
                'ListSize',[420 400],'ListString', str) ; 
 
epilogh = d(s).name ;            
% load stereo sound
[hxhtiko,Fs] = audioread(epilogh) ;

cd ..

% stereo2mono
if( size( hxhtiko, 2)>1)
    hxhtiko = (hxhtiko(:,1) + hxhtiko(:,2))./2 ;
    %hxhtiko2 = (hxhtiko2(:,1) + hxhtiko2(:,2))./2 ;
end

% normalise input
hxhtikoN = normalisedVEC(hxhtiko) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------- Feed-Forward Side Chain ---------------------------

% abs-value 
hxhtikoLD = max( abs(hxhtikoN), 1e-6) ;

% linear to dB conversion
hxhtiko_dB = 20*log10(hxhtikoLD) ;

% choice -> apanthsh tou xrhsth meta apo orismo idiothtwn compressora
choice = '' ;
while( ~strcmp(choice,'Yes') )

    prompt = {'Threshold(dB):','Ratio:','Attack Time(sec.):',...
        'Release Time(sec.):','Knee Width(dB):','Level Detection(rms/peak):'} ;
    dlg_title = 'Compressor ' ;
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

    % Threshold
    %Threshold = -45 ;

    % Ratio
    %Ratio = 8 ;

    % Knee Width
    %KneeW = 10 ;

    inputEX = -60:5:0 ;
    inputEX = inputEX' ; 
    
    if(isequal(KneeW,0))
        outputEXSG(inputEX<=Threshold) = inputEX(inputEX<=Threshold) ;
        outputEXSG(inputEX>Threshold) = Threshold + (inputEX(inputEX>Threshold) - Threshold)/Ratio ;
    
    else
    
        for(i=1:length(inputEX))
    
            % case1 : 2*(xG - T) < - W
            if( 2*(inputEX(i,1)-Threshold) < -KneeW )
                outputEXSG(i,1) = inputEX(i,1) ;
            % case2: 2*abs(xG-T) <= W    
            elseif( 2*abs(inputEX(i,1)-Threshold) <= KneeW )
                outputEXSG(i,1) = inputEX(i,1)+(((1/Ratio)-1)*((inputEX(i,1)-Threshold+KneeW/2)^2))/(2*KneeW) ;
            % case3: 2*(xG - T) > -W
            elseif( 2*(inputEX(i,1)-Threshold) > KneeW )
                outputEXSG(i,1) = Threshold + (inputEX(i,1) - Threshold)/Ratio ;    
            end
    
        end

     end
        
     figure(1)
     plot( inputEX, outputEXSG)
     hold on ;
     plot( inputEX, inputEX)
     plot( Threshold*ones( size(inputEX,2), size(inputEX,1)),inputEX)
     hold off ;
     xlabel('Input (dB)') ;
     ylabel('Output (dB)') ;
     legend(['Compression with Ratio=',num2str(Ratio),':1'],'No Compression(1:1)','Threshold') ;
     title('Compressor Static Characteristic Curve') ;
     
     options.Interpreter = 'tex';
     % Include the desired Default answer
     options.Default = 'Yes';
     % Use the TeX interpreter in the question
     qstring = 'Do you want to apply the Compressor with these settings ?';
     choice = questdlg( qstring,'Apply Compressor',...
          'Yes','Change Settings', options) ;
    
end

% attack time coefficient
Co_att = exp(-log(9)/(Fs*T_att)) ;
%Co_att = exp((-log(9))/(Fs*T_att)) ;
% release time coefficient
Co_rel = exp(-log(9)/(Fs*T_rel)) ;
%Co_rel = exp((-log(9))/(Fs*T_rel)) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------- Feed-Forward Side Chain ---------------------------


% linear to dB conversion
hxhtiko_dB = 20*log10(hxhtikoLD) ;

inputG = hxhtiko_dB ; 

outputG = zeros( size(inputG) ) ;
coga = outputG ;
% hard knee characteristic - hard half-wave rectifier
if(isequal(KneeW,0))
    
    for(i=1:length(inputG'))
    
        if(inputG(i,1)<=Threshold)
        
            outputG(i,1) = inputG(i,1) ;
        else
            
            outputG(i,1) = Threshold + (inputG(i,1) - Threshold)/Ratio ;      
        end
             
        % computer gain on soft knee characteristic
        coga(i,1) = -outputG(i,1) + inputG(i,1) ;
    end
    
else
            
     for(i=1:length(inputG'))
     % soft knee characteristic - soft half-wave rectifier
    
        % case1 : 2*(xG - T) < - W
        if( 2*(inputG(i,1)-Threshold) < -KneeW )
            outputG(i,1) = inputG(i,1) ;
        % case2: 2*|(xG-T)| <= W    
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


% -------------------- level detection-gain smoothing ---------------------

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


% -------------------- Make-Up Gain Addition ------------------------------

% auto make up gain computation - negative of the computed gain for 0dB

if((KneeW/2)<Threshold)
    tempMUG = 0 ;
    %c1 = c1 + 1 ;
elseif((2*abs(Threshold))<=KneeW)
    tempMUG = ((1/Ratio - 1)*(Threshold-(KneeW/2))^2)/(2*KneeW) ;
    %c2 = c2 + 1 ;
elseif(Threshold<(-KneeW/2))
    tempMUG = Threshold + Threshold/Ratio ;
    %c3 = c3 + 1 ;
end

magain = -smga - tempMUG ;


% --------------------------- dB 2 Linear --------------------------------- 

% the gain reduction to be applied to the signal
lingain = 10.^(magain./20) ;
ThresLin = 10^(Threshold/20) ;

% ------------------------ Compressed Output ------------------------------
compOUT = hxhtikoN.*lingain ;

% normalised output
compOUTn = normalisedVEC(compOUT) ;

figure(2)
subplot(1,2,1)
plot(hxhtikoN)
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

figure(3)
plot(hxhtikoN)
hold on ;
plot(compOUTn)
plot(ThresLin*ones(size(hxhtiko)))
hold off ;
legend('Input','Output','Threshold') ;
title('Input/Output Plot') ;
xlabel('Samples(n)') ;
ylabel('Amplitude') ;

% % emfanish parathyrou sto xrhsth gia akroash apotelesmatos tou sympiesth 
while(1)

    options.Interpreter = 'tex';
    % Include the desired Default answer
    options.Default = 'After';
    % Use the TeX interpreter in the question
    qstring = 'Do you want to listen to the audio sample before or/and after the Compressor appliannce ?';
    choiceC = questdlg( qstring,'Compressor Audio Result',...
        'Before','After',  'Again', options) ;
    
    if(strcmp(choiceC,'Before'))
       sound( hxhtiko, Fs) ;
    
    elseif(strcmp(choiceC,'After'))
        sound( compOUTn, Fs) ;
    
    else
        close all ; 
        break ;
    
    end

end

%choiceC = '' ;

end
