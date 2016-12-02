function [cand_pss_ca, cand_pss_ca_order] = Pss_reliability_filtering(P2xcor_num22,Ppss_th)

global cs_para;


cand_pss_ca = []; 
cand_pss_calc_order = [];

for i = 1:3
    tmpindex = find(P2xcor_num22(i,:) >= Ppss_th(i,1));
    if (size(tmpindex,2) ~= 0)
        cand_pss_calc(:,2) = tmpindex;
        cand_pss_calc(:,3) = P2xcor_num22(i,cand_pss_calc(:,2));
        cand_pss_calc(:,4) = 0; % default
        cand_pss_calc(1:size(cand_pss_calc,1),1) = i;
        
        % 在cand_pss_calc中计算各个pss的可靠峰
%         cand_pss_calc_order = [1 2954 202; 1 2960 80; 1 2983 40];  % Unit test
        cand_pss_calc_order = flipud(sortrows(cand_pss_calc,3));
        if (cs_para.pss_reliability_filter_flag == 1)
            WL_pss_peak_det = 5;      % 主峰左右各多少点纳入监测集
            outwindow_index = find((cand_pss_calc_order(:,2) > cand_pss_calc_order(1,2)+WL_pss_peak_det) | (cand_pss_calc_order(:,2) < cand_pss_calc_order(1,2)-WL_pss_peak_det));
            if (size(outwindow_index,1) == 0)||...
                    ((cand_pss_calc_order(1,3)/ cand_pss_calc_order(outwindow_index(1,1),3)>= cs_para.pss_reliability_1st_out1st_ratio) &&...
                    (cand_pss_calc_order(1,3)/Ppss_th(i,1) >= cs_para.pss_reliability_1st_thresold_ratio))
                cand_pss_calc_order(1,4) = 1;
                cand_pss_calc_order(outwindow_index,:) = [];
            end;
        end;
        
        cand_pss_ca = [cand_pss_ca;cand_pss_calc_order];
        cand_pss_calc = [];
    else
        continue;
    end;   % end
end;         
cand_pss_ca_order = flipud(sortrows(cand_pss_ca,3));
clear cand_pss_calc cand_pss_calc_order i; 


