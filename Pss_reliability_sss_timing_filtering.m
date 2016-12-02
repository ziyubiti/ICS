function [cand_pss_ca_order2] = Pss_reliability_sss_timing_filtering(cand_pss_ca_order,TimeDataDwsampling,SF_num)

global cs_para;

cand_pss_ca_order2 = cand_pss_ca_order;

if (cs_para.pss_reliability_sss_timing_filter_flag == 1)
    cand_pss_ca_order2(:,5) = 0;
    for j = 1:size(cand_pss_ca_order2,1)
        NID2No = cand_pss_ca_order2(j,1);
        NID2 = NID2No-1;
        
        %% SSS timing for pss reliability filter
        
        Pos_pss_timing = cand_pss_ca_order2(j,2); %pos(NID2No,1);  %randi(4800,1,1);%
        Pos_sss_timing_ori = Pos_pss_timing - 4544/32 - 2048/32;
        if (Pos_sss_timing_ori<0)
            Pos_sss_timing_ori = Pos_sss_timing_ori + 4800;
        end;
        WL_sss_timing = 64;
        PQsss_sum = 0;
        PQsss_sum2 = 0;
        for iss = 1:SF_num/5
            Pos_sss_timing = Pos_sss_timing_ori + (iss-1)*4800;
            Pos_sss_timing_start = Pos_sss_timing - WL_sss_timing;
            Pos_sss_timing_end = Pos_sss_timing + WL_sss_timing;
            Tsss_timing = TimeDataDwsampling(Pos_sss_timing_start:Pos_sss_timing_end + cs_para.N_FFT-1,1);
            
%             t = load('D:\00 work\10 算法仿真\90 数值处理\time_data0102_dec.txt');
% t2 = t(1:2:end,1)+1j*t(2:2:end,1);
% t3 = t2/2.^14;
%         Tsss_timing = t3;
        
            for sss_timing_index = 1:2*WL_sss_timing
                tmp_reserve = flipud(Tsss_timing([sss_timing_index+33:sss_timing_index+cs_para.N_FFT-1],1));
                Qsss(sss_timing_index,1) = sum(Tsss_timing(sss_timing_index+1:sss_timing_index+31,1).*tmp_reserve);
            end;
            PQsss = abs(Qsss).*abs(Qsss);
            PQsss_sum = PQsss_sum + PQsss;
            
            if (cs_para.sss_timing_slide_flag == 1)
                WLslide = cs_para.sss_timing_slide_WL;
                PQsss_tmp = [ones(1,WLslide/2)*1e-10  PQsss.'  ones(1,WLslide/2)*1e-10];
                for jsss = (WLslide/2)+1:size(PQsss_tmp,2)-(WLslide/2)
                    PQsss_tmp2(1,jsss) = PQsss_tmp(1,jsss)/mean(PQsss_tmp(1,[jsss-WLslide/2:jsss-1 jsss+1:jsss+WLslide/2]),2);
                    %                 if (xcortmp2(i,jpss) == NaN)
                    %                     xcortmp2(i,jpss) = 0;
                    %                 end;
                end;
                PQsss_tmp3(1,:) = PQsss_tmp2(1,[(WLslide/2)+1:size(PQsss_tmp,2)-(WLslide/2)]);
                PQsss_sum2 = PQsss_sum2 + PQsss_tmp3.';
            end;
            
            
        end;  % end of sf loop
        
        if (cs_para.sss_timing_slide_flag == 1)
            PQsss_sum3 = PQsss_sum2;
        else
            PQsss_sum3 = PQsss_sum;
        end;
        noise_power_sss_timing = mean(PQsss_sum3);
        Psss_timing_th = cs_para.sss_timing_threshold*noise_power_sss_timing;
        
        [Psss_timing_peak Psss_timing_peak_index] = max(PQsss_sum3);
        
        if (Psss_timing_peak >= Psss_timing_th) && ((Psss_timing_peak_index == 64) || (Psss_timing_peak_index == 65) || (Psss_timing_peak_index == 66))
            cand_pss_ca_order2(j,5) = 1;
        end;
        
        figure();
        plot([0 size(PQsss_sum3,1)],[Psss_timing_th Psss_timing_th],'r',[0 size(PQsss,1)],[noise_power_sss_timing  noise_power_sss_timing],'c',[1:size(PQsss_sum3,1)],PQsss_sum3,'b');
        grid on;
%         title(int2str(DelayT));
        set (gcf,'Position',[1,500,400,300], 'color','w')
        
        %   end;  % end of vary Delay, Unit Test
        
        
    end;   % end of cand_pss_ca_order loop for sss timing filter
    
    % update the cand_pss_ca_order, remove the pss  not satisfied sss timing
    tmpindex = find(cand_pss_ca_order2(:,5) == 0);
    cand_pss_ca_order2(tmpindex,:) = [];
 
    
end;  % end of sss timing filter

