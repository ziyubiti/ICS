             


function [phy_sc_index] = phy_scmap(N_FFT,N_RB,Nsc_RB,MA_option)
switch MA_option
    case {'DFT-S-OFDMA'}
        phy_sc_index=[(N_FFT/2-N_RB*Nsc_RB/2+1):(N_FFT/2+N_RB*Nsc_RB/2)];
    case {'OFDMA'}
        phy_sc_index=[(N_FFT/2-N_RB*Nsc_RB/2+1):N_FFT/2,(N_FFT/2+2):(N_FFT/2+N_RB*Nsc_RB/2+1)];
    otherwise
     error('multiple access option not defined at present.');
end