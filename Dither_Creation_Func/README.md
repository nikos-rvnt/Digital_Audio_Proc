Function Generating Dither Noise to be added to a signal. 

                                y = x + d + v ,
                                
where x -> inout signal, v -> quantisation noise and d -> dither noise that gives a better noise to quant noise.
                                
Dither can be used as added noise in order to cancel bad effects of quantisation or after tha mixing and during the master process
to give better audible to human result.

In this function are implemented 4 different pdf-shape of dither, for different dither input flag.

For flag == 1 -> Rectangular pdf,
            2 -> Triangular pdf,
            3 -> High-Pass triangular pdf,
            4 -> Gaussian pdf,
            
E.g. calling:

              dither = DitherGen(DitherFlag,LSB,npoints, InSignal) ;
              
returns a dither noise signal of equal length with InSignal (which is x). LSB is calculated as :

LSB = (up_quant_limit - down_quant_limit)/2^N-1 , where N is number of bit resolution. npoints -> length of InSignal.
