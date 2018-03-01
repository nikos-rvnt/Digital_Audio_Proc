function [dither] = DitherGen( DitherFlag, LSB, npoints, xsignal)

    % ---------------------------------------------------------------------
    % DitherGen() - Dither Noise Generation Function, where :
    %
    % Input :
    %    - DitherFlag : 1 -> RPDF
    %                   2 -> TPDF
    %                   3 -> High-Pass TPDF
    %                   4 -> GaussianPDF
    %
    %    - LSB : PCM quantisation step
    %
    %    - npoints : length of quantisation noise signal
    %
    %    - xsignal : input signal to be quantised
    %
    % Output :
    %    - dither : the created dither signal to be added to the
    %               quantisation noise and the input signal.
    % ---------------------------------------------------------------------
    % x - input signal, v - quantisation noise, d - dither signal and x, d, 
    % v independent of each other, then :
    %
    %                       y = x + d + v ;
    %
    % where y - quantized input signal.
    % ---------------------------------------------------------------------

    % generate random noise (white noise) in range of 0 to 1(exclusive)
    noise = rand( size( xsignal, 1), size( xsignal, 2)) ;
    
    % ---------------- RPDF - rectangular pdf dither noise ----------------
    if isequal( DitherFlag, 1)
        % uniformly distributed & independent of the input signal 
        noise = -LSB/2 + LSB.*noise ;
    
    % --------------------- TPDF - triangular -----------------------------   
    elseif isequal( DitherFlag, 2)
        % addition of 2 uniformly distributed & independent of each other
        % and of the input signal
        noise = -LSB + 2*LSB.*(noise + rand(size(noise))) ; 
    
    % --------------------- High - Pass TPDF ------------------------------
    elseif isequal( DitherFlag, 3)
        % substraction of every sample value of uniformly distributed noise
        % by its previous sample value
        noise_pre = 0 ;
        noiseH = zeros( size(noise) ) ;
        for i=1:length(noise)
           % dt(i) = du(i) - du(i-1) 
           noiseH(i) = noise(i) - noise_pre ;
           noise_pre = noise(i) ;
        end
        noise = - LSB + 2*LSB.*noiseH ;
        
    % --------------------- GPDF - gaussian -------------------------------   
    elseif isequal( DitherFlag, 4)
        % independent white gaussian noise with power = Ïƒ >= A/2
        noise = wgn( size(noise), 20*log10((LSB^2)/3)) ;         
    
    % unsupported dither type
    else 
        error('Unsupported dither type! Choose among values: 1,2 or 3.') ;
        
    end
    
    dither = noise ;

end
