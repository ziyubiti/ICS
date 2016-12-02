
global sys_para;

sys_para.N_RB_DL = 100;
sys_para.Nsc_RB = 12;
sys_para.N_FFT = 2048;
sys_para.N_DL_symb = 7;
sys_para.sub_carrier_spacing = 15000;
sys_para.fs = sys_para.sub_carrier_spacing*sys_para.N_FFT;
sys_para.ncp_length = [160 144];

