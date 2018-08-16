## Dynamic Range Compression

Project implementation for *Digital Audio Technology* class of MSc programm [*SIGNAL PROCESSING AND COMMUNICATION SYSTEMS-SPCOMS*](http://xanthippi.ceid.upatras.gr/dsp/). 

* Implemented on Matlab R2015b. 
* For the digital range compressor algorithm implementation, I was based on the next paper, ["Digital Range Compressor Design - A Tutorial and Analysis" D. Giannoulis, M. Massberg, J. Reiss](https://pdfs.semanticscholar.org/f1b2/0a5681e6ef7080e5b5fbce81911c6873543c.pdf).

As seen below, the process of range compression is described by a **feed forward side-chain block with an auto make-up gain addition** before output. 

![ ](https://raw.githubusercontent.com/nikos-rvnt/Digital_Audio_Processing/master/Dynamic_Range_Compression/drc_algo_block_diagram.jpg)

The basic matlab code is *"Digital_Range_Compression.m"*, which while running it some user-input windows are poped up, in order for a music sample to be chosen (about 6 sec. cropped .wav files), or for setting up the compressor characteristics. The program ends if the user decides so. Hence, when the compression process finishes, some plots come up comparing tha music sample before and after compression and a window asking the user if he want to listen to the sample before or/and after the compression and lastly if the user wants to repeat the process for another sample or for the same with different compressor parameters. 

The above described process can be better understood from the pictures below:

  - Choosing the audio sample:
  
    ![ ](https://raw.githubusercontent.com/nikos-rvnt/Digital_Audio_Processing/master/Dynamic_Range_Compression/comp_1.png)

  - Setting the compressor characteristics:
  
    ![ ](https://raw.githubusercontent.com/nikos-rvnt/Digital_Audio_Processing/master/Dynamic_Range_Compression/comp_2.png)
  
  - The compressor curve and a question to the user to proceed with this settings:
  
    ![ ](https://raw.githubusercontent.com/nikos-rvnt/Digital_Audio_Processing/master/Dynamic_Range_Compression/comp_3.png)

  - The audio sample before and after the compressor process:
  
    ![ ](https://raw.githubusercontent.com/nikos-rvnt/Digital_Audio_Processing/master/Dynamic_Range_Compression/comp_4.png)

    ![ ](https://raw.githubusercontent.com/nikos-rvnt/Digital_Audio_Processing/master/Dynamic_Range_Compression/comp_5.png)


  - A question to the user whether he wants to listen to the compressed audio sample, or he wants to repeat the process with the same or another sample (pressing X the program ends):

    ![ ](https://raw.githubusercontent.com/nikos-rvnt/Digital_Audio_Processing/master/Dynamic_Range_Compression/comp_6.png)


On the other hand, *"compTest.m"* file and *"DRC_test.m"* function respectively, implement a test on the functionality of the compressor giving as input a vector of different sinusoidals one after the other sampled on **44100 Hz** each one with fund. frequency **4000 Hz** and amplitudes **.24, .92, .48** respectively. 

Lastly, *"DRC_figures.m"* is just a code showing the change on the signal before and after each stage of the feed-forward side chain block compressor. 
