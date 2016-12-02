

function [freq_domain_signal]=ofdm_fft(time_domain_signal,N_FFT,dim)

% sqrt(N_FFT) Scale the signal power to constant for MATLAB fft function 
freq_domain_signal=1/sqrt(N_FFT)*fft(time_domain_signal,N_FFT,dim);
% Rearrange the data order for fft 
freq_domain_signal=fftshift(freq_domain_signal,dim);

% return;