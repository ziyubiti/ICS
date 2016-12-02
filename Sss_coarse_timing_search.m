% delayT仅用于单元测试时画图时使用，不影响主程序功能。

function [Peak_sss_timing_value,  Peak_sss_timing_idx,noise_power_sss_coarse_timing] = Sss_coarse_timing_search(TimeDataDwsampling,DelayT)

global cs_para;


Tsss_timing_cand = [TimeDataDwsampling;TimeDataDwsampling(1:64,1)];

for sss_timing_index_all = 1:size(TimeDataDwsampling,1)
    tmp_reserve = flipud(Tsss_timing_cand([sss_timing_index_all+33:sss_timing_index_all+cs_para.N_FFT-1],1));
    Qsss_all(sss_timing_index_all,1) = sum(Tsss_timing_cand(sss_timing_index_all+1:sss_timing_index_all+31,1).*tmp_reserve);
end;
PQsss_all = abs(Qsss_all).*abs(Qsss_all);


if (cs_para.sss_timing_slide_flag == 1)
    WLslide = cs_para.sss_timing_slide_WL;
    PQsss_tmp = [ones(1,WLslide/2)*1e-10  PQsss_all.'  ones(1,WLslide/2)*1e-10];
    for jsss = (WLslide/2)+1:size(PQsss_tmp,2)-(WLslide/2)
        PQsss_tmp2(1,jsss) = PQsss_tmp(1,jsss)/mean(PQsss_tmp(1,[jsss-WLslide/2:jsss-1 jsss+1:jsss+WLslide/2]),2);
        %                 if (xcortmp2(i,jpss) == NaN)
        %                     xcortmp2(i,jpss) = 0;
        %                 end;
    end;
    PQsss_tmp3(1,:) = PQsss_tmp2(1,[(WLslide/2)+1:size(PQsss_tmp,2)-(WLslide/2)]);
    PQsss_all = PQsss_tmp3.';
end;
PQsss_shape = reshape(PQsss_all,4800,[]);
PQsss_sum = sum(PQsss_shape,2);

[Peak_sss_timing_value  Peak_sss_timing_idx] = max(PQsss_sum);
noise_power_sss_coarse_timing = mean(PQsss_sum);
Psss_coarse_timing_th = cs_para.sss_timing_threshold*noise_power_sss_coarse_timing;


figure();
plot([0 size(PQsss_sum,1)],[Psss_coarse_timing_th Psss_coarse_timing_th],'r',[0 size(PQsss_sum,1)],[noise_power_sss_coarse_timing  noise_power_sss_coarse_timing],'c',[1:size(PQsss_sum,1)],PQsss_sum,'b');
grid on;
title(strcat('SSS timing,DelayT = ',int2str(DelayT)));
set (gcf,'Position',[1,100,400,300], 'color','w')
clear Tsss_timing_cand sss_timing_index_all tmp_reserve Qsss_all PQsss_tmp jsss PQsss_tmp2 PQsss_tmp3 PQsss PQsss_shape PQsss_sum;

