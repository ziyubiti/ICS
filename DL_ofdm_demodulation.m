% filename:             ofdm_Demodulation_downlink
% Description:          Remove CP, FFT operation
%
%
% Input:                receive_frame_signals             2D :Rx*30720
%                       receive_frame_signals0            
% Output:               receive_antenna_symbols           3D :rxa*1200*14
%                       receive_antenna_symbols0          
% Calling:

%
function [receive_antenna_TD receive_antenna_TD0 receive_antenna_symbols,receive_antenna_symbols0,UT_FD_all] =  DL_ofdm_demodulation(receive_frame_signals,receive_frame_signals0)

global sys_para;


[rxa_num,SAMPLES_NUM] = size(receive_frame_signals);
N_RB_DL = sys_para.N_RB_DL;
Nsc_RB = sys_para.Nsc_RB;
N_FFT = sys_para.N_FFT;
N_DL_symb = sys_para.N_DL_symb;
% sc_spacing=sys_para.sub_carrier_spacing;
normal_cp_length = sys_para.ncp_length;
% samples_num=N_FFT*sc_spacing/1000;%表示一个子帧内的抽样点数
N_SYMBLS = N_DL_symb*2;
phy_sc_index = phy_scmap(N_FFT,N_RB_DL,Nsc_RB,'OFDMA');
N_CP1 = normal_cp_length(1);
N_CP2 = normal_cp_length(2);


% 处理一个时隙. 上下行的时隙结构不同，应该分离开来
for n_rxa=1:1:rxa_num
    subframe_signals_tem  = squeeze(receive_frame_signals(n_rxa,:));
    subframe_signals0_tem = squeeze(receive_frame_signals0(n_rxa,:));
    data_dim = 1;
    % 去掉循环前缀
    symbol_leng = 0;
    for n_symbols = 1:N_SYMBLS
        if (n_symbols == 1)|| (n_symbols==8)
            symbol_leng_tem1 = N_FFT+N_CP1;
            receive_frame_signals_tem = subframe_signals_tem(symbol_leng+1:symbol_leng+symbol_leng_tem1);
            receive_frame_signals_tem0 = subframe_signals0_tem(symbol_leng+1:symbol_leng+symbol_leng_tem1);
            subframe_signals_nocp(:,n_symbols)  = cp_remove(receive_frame_signals_tem,N_CP1,data_dim);
            subframe_signals_nocp0(:,n_symbols)  = cp_remove(receive_frame_signals_tem0,N_CP1,data_dim);
            symbol_leng = symbol_leng + symbol_leng_tem1;
        else
            symbol_leng_tem2 = N_FFT+N_CP2;
            receive_frame_signals2_tem = subframe_signals_tem(symbol_leng+1:symbol_leng+symbol_leng_tem2);
            receive_frame_signals2_tem0 = subframe_signals0_tem(symbol_leng+1:symbol_leng+symbol_leng_tem2);
            subframe_signals_nocp(:,n_symbols)  = cp_remove(receive_frame_signals2_tem,N_CP2,data_dim);
            subframe_signals_nocp0(:,n_symbols)  = cp_remove(receive_frame_signals2_tem0,N_CP2,data_dim);
            symbol_leng = symbol_leng + symbol_leng_tem2;
        end
    end
    % 变换到频域
    subframe_symbols   = ofdm_fft(subframe_signals_nocp,N_FFT,data_dim);
    subframe_symbols0  = ofdm_fft(subframe_signals_nocp0,N_FFT,data_dim);
    
    UT_FD_all(n_rxa,:,:) = subframe_symbols;
    %
    subframe_symbols   = subframe_symbols([phy_sc_index],:);%有用的子载波
    subframe_symbols0  = subframe_symbols0([phy_sc_index],:);

    receive_antenna_TD(n_rxa,:,:) = subframe_signals_nocp;
    receive_antenna_TD0(n_rxa,:,:) = subframe_signals_nocp0;
    receive_antenna_symbols(n_rxa,:,:) = subframe_symbols;
    receive_antenna_symbols0(n_rxa,:,:) = subframe_symbols0;
end % one antenna over
