% 由于降采样后数据较少；利用原时域数据的cp最后72个点进行FFO估计；
% 根据估计的FFO对降采样后数据进行补偿。



function [TimeDataDwsampling_ffoCP, foe_cp_ffo] = Foe_cp_ffo_estimation(TimeData,TimeDataDwsampling,Peak_sss_timing_idx,SF_num)

global cs_para;
global sys_para;


Peak_sss_timing_idx_up_orig = (Peak_sss_timing_idx-1)*cs_para.Ratio + 1;
xcor_sum = 0;
for i = 1:SF_num/5
    Peak_sss_timing_idx_up = Peak_sss_timing_idx_up_orig + 153600*(i-1);
    data1 = TimeData(Peak_sss_timing_idx_up-72:Peak_sss_timing_idx_up-1,1);
    data2 = TimeData(Peak_sss_timing_idx_up-72+sys_para.N_FFT:Peak_sss_timing_idx_up-1+sys_para.N_FFT,1);
    xcor_sum = xcor_sum + sum(data2.*conj(data1));
end

foe_cp_ffo = angle(xcor_sum)/2/pi*15000;
comp_phase0 = exp(-j*2*pi*foe_cp_ffo./(cs_para.N_FFT*sys_para.sub_carrier_spacing).*(0:size(TimeDataDwsampling,1)-1));
TimeDataDwsampling_ffoCP = TimeDataDwsampling .*comp_phase0.';


