%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                   DPMS SESE  -  Psifiakh Texnologia Hxou                %
%                                                                         %
%           Ylopoihsh Sympiesth Dynamikhs Perioxhs (Compressor)           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

fs =   44100 ;
f_sin = 4000 ;

t_sin = 3 ;
t = 1:(1/fs):(t_sin-1/fs) ;


% 3 sinusoidal functions generation
sin1 = 0.24*sin( 2*pi*t*f_sin ) ;
sin2 = 0.92*sin( 2*pi*t*f_sin ) ;
sin3 = 0.48*sin( 2*pi*t*f_sin ) ;

% 3 sinusoidals v 
sinus = [ sin1, sin2, sin3] ;

% sinus = audioread('input_ampSweep_-56_-6dB.wav') ;
% sinus = sinus' ;

% dhmiourgia parathyrou hanning
winHam = hamming( max(size(sinus,1),size(sinus,2)) ) ;
% parathyrwsh hmitonou
sinusW = sinus.*winHam' ;

% compressor parameters
Threshold = -25 ;
KneeW =  5 ;
Ratio = 10 ;
T_att = 0.0015 ; % 0.0083 ;
T_rel = 0.0180 ; % 0.0074 ;
LVLD = 'peak' ;

% kanonikopoihsh sto [-1,1]
%sinusW  = normalisedVEC(sinusW) ;
sinus  = normalisedVEC(sinus) ;

% Digital Range Compressor Analysis
[ sinusCOMP, smgaS, coga] = DRC_test( sinus, Threshold, KneeW, Ratio, T_att, T_rel,LVLD,fs) ;

% dhmiourgia parathyrou hanning
winHam = hamming( max(size(sinusCOMP,1),size(sinusCOMP,2)) ) ;
% parathyrwsh hmitonou
sinusCOMP = sinusCOMP.*winHam' ;

% kanonikopoihsh sto [-1,1]
%sinusWN = normalisedVEC(sinusW) ;
sinusCOMPN = normalisedVEC(sinusCOMP) ;

% grammiko katwfli
TLin = 10^(Threshold/20) ;

figure(1)
subplot(2, 2, 1)
plot(sinus)
hold on ;
plot( TLin*ones(size(sinus)))
hold off ;
xlabel('samples(n)') ;
ylabel('Amplitude') ;
title('Input Signal') ;
subplot(2,2,2)
plot(coga)
xlabel('samples(n)') ;
ylabel('Amplitude(dB)') ;
title('Computed Gain') ;
subplot(2,2,3)
plot(smgaS)
xlabel('samples(n)') ;
ylabel('Amplitude(dB)') ;
title('Smoothed Gain') ;
subplot(2,2,4)
plot(sinusCOMPN)
hold on ;
plot( TLin*ones(size(sinus)))
hold off ;
xlabel('samples(n)') ;
ylabel('Amplitude') ;
title('Compressed Signal') ;

figure(2)
subplot(3,1,1)
plot( sinus )
hold on ;
plot( sinusCOMPN )
plot( TLin*ones(size(sinusW)))
hold off ;
xlabel('samples') ;
ylabel('Amplitude') ;
legend('Input','Output','Threshold') ;
title('Before & After Compressor Appliance') ;

subplot(3,1,2)
plot( sinus )
hold on ;
plot( TLin*ones(size(sinusW)) )
hold off ;
xlabel('samples') ;
ylabel('Amplitude') ;
title('Sine Sequence - Compressor Input') ;
subplot(3,1,3)
plot( sinusCOMPN )
hold on ;
plot( TLin*ones(size(sinusW)) )
hold off ;
xlabel('samples') ;
ylabel('Amplitude') ;
title('Sine Sequence - Compressor Output') ;

% ---------------------- Ypologismos THD & SNR ----------------------------

% ypologismos fasmatos gia eisodo sinusW
spec1 = fft( sinusW ) ; % spectrum
magSpec1 = abs( spec1 ); % spectral magnitude
magSpec_dB1 = 20*log10(magSpec1/max(abs(spec1))); % spectral magnitude in dB-FS
freq1 = 1:44100 ;
% ypologismos thesis themeliwdous syxnothtas
fc1 = find(magSpec1==max(magSpec1)) ;
% frequency vector
l = length(sinusW);    
fvec = 0:fs/l:fs-1/l; 
% ypologismos armonikhs paramorfwshs+thoryvou se dB
[ THD_NOISE ] = 10*log10( CalcTHDN( magSpec1, fvec(1:l/2), fc1(1,1))) ;

% ypologismos fasmatos gia eksodo sinusCOMP
spec2 = fft( sinusCOMP ) ; % spectrum
magSpec2 = abs( spec2 ) ; % spectral magnitude
magSpec_dB2 = 20*log10(magSpec2/max(abs(spec2))); % spectral magnitude in dB-FS
freq2 = 1:44100 ;
% ypologismos thesis themeliwdous syxnothtas
fc2 = find(magSpec2==max(magSpec2)) ;
% frequency vector
l = length(sinusW);    
fvec = 1:fs/l:fs-1/l; 
% ypologismos armonikhs paramorfwshs+thoryvou se dB
[ THD_NOISE2 ] = 10*log10( CalcTHDN( magSpec2, fvec(1:l/2), fc2(1,1))) ;

% -------------------- Elegxos Pistothtas Perivallousas -------------------

% ypologismos perivallousas me filtro Hilbert
inpEnv = abs( hilbert(sinus) ) ;
outEnv = abs( hilbert(sinusCOMPN) ) ;

% ypologismos dianysmatos eterosysxetishs eisodou-eksodou sympiesth
IOEnv = corr( inpEnv', outEnv') ;

figure(3)
plot(inpEnv)
hold on ;
plot(outEnv)
plot( TLin*ones(size(sinus)) )
hold off ;
legend('Input Envelope','Output Envelope','Threshold') ;
xlabel('samples(n)') ;
ylabel('Abs Amplitude') ;
title(['Input & Output Envelope with Correlation:',num2str(IOEnv)]) ;
