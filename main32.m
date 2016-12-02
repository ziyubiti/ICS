
clear;
close all;
% 32倍降采样处理，支持DSP 65阶滤波器设置。


global cs_para;
global sys_para;

LNAGain = 8.0;
AD9361Gain = 28;
TestVectorPrintFlag = 0;

sys_para_set;
cs_para_set;
SF_num = 40;                                   % 用于评估计算的子帧数；  

DelayT = 0; % mod(randi(33,1,1),33)-16;                                               % 采样点偏差，-DownSamplingRatio到DownSamplingRatio, 16倍下采样时，为-16到+16
% 设为正值表示延时采样信号，延时后信号相关的峰值位置提前；设为负值表示提前采样信号，这时峰值延后。

%% load single cell Test Vector IQ data
     for DelayT = 0%0:4:16;       % unit test
%    for deltaF = -7.5e3:1e3:7.5e3          % unit test
% %  Tdata = load('D:\00 work\03 测试\02 系统测试\2013\ZCTT_tx_signal_time_tdd1_10sf.txt');
%   Tdata = load('20140704ICS\cellID3（0 1).txt');%; 
%    Tdata = load('20141229FDD\Case 4\ICS_40ms_mutli_FDD_Case4A.dat');
%   D = Tdata(:,1)+1j*Tdata(:,2);  clear Tdata;
% DelayT = 0;%mod(randi(33,1,1),33)-16;  ;                                            % 设为正值表示延时采样信号，延时后信号相关的峰值位置提前；设为负值表示提前采样信号，这时峰值延后。
% deltaF = 0e3;
% n = 1:30720*SF_num;
% foe = exp(j*2*pi*deltaF./(2048*15000).*(n-1));
% D = D.*foe.'+ eps;  
%  TimeData40sf = D + eps;

%% load multi cell test vector
% Data_cell1 = load('20140704ICS\cellID3（0 1).txt'); 
% Data_cell2 = load('20140704ICS\cellID106(1 35).txt'); 
% Data_cell3 = load('20140704ICS\cellID344(2 114).txt'); 
% Data_cell4 = load('20140704ICS\cellID423(0 141).txt'); 
% D(:,1) = Data_cell1(:,1) + 1j*Data_cell1(:,2);          
% D(:,2) = Data_cell2(:,1) + 1j*Data_cell2(:,2); 
% D(:,3) = Data_cell3(:,1) + 1j*Data_cell3(:,2); 
% D(:,4) = Data_cell4(:,1) + 1j*Data_cell4(:,2); 
% PowerOffset = sqrt(10.^([0 -5 -10 -15]/10));
% deltaF = [0 ; -2.6e3; 3.5e3; 7.8e3];
% % deltaF = [0 ; (mod(randi(16,1,1),16)-7.5)*1e3; (mod(randi(16,1,1),16)-7.5)*1e3; (mod(randi(16,1,1),16)-7.5)*1e3];
% TimeOffset = [0 10 30 90];
% % TimeOffset = [0  randi(160,1,1)  randi(160,1,1)  randi(160,1,1)];
% n = 1:307200;
% foe = exp(j*2*pi*deltaF./(2048*15000)*(n-1));
% D = D.*repmat(PowerOffset,307200,1).*foe.';
% D2 = zeros(max(TimeOffset)+size(D,1),4);
% D2(max(TimeOffset)+1:size(D,1)+max(TimeOffset),:) = D;
% Data_Uu = D2(max(TimeOffset)-TimeOffset(1)+1:max(TimeOffset)-TimeOffset(1)+307200,1)...
%           + D2(max(TimeOffset)-TimeOffset(2)+1:max(TimeOffset)-TimeOffset(2)+307200,2)...
%           + D2(max(TimeOffset)-TimeOffset(3)+1:max(TimeOffset)-TimeOffset(3)+307200,3)...
%           + D2(max(TimeOffset)-TimeOffset(4)+1:max(TimeOffset)-TimeOffset(4)+307200,4);
% clear Data_cell1 Data_cell2 Data_cell3 Data_cell4 D D2 foe n ;
% TimeData40sf = Data_Uu + eps;


if(TestVectorPrintFlag)
    fid3 = fopen('ICS_10ms_pci3_6.8k.dat','w+');
    for i = 1:1:size(TimeData40sf)
        fprintf(fid3,'%10.6f     %10.6f \n' ,real(TimeData40sf(i,1)),imag(TimeData40sf(i,1)));
    end;
    
    fclose('all');
end;


%% load DSP data   
Tdata = load('D:\2016 ver\0826 cqi\IQ_data\TD_LowHigh_40ms.txt');
TimeData40sf = (Tdata(1:2:end,1)+1j*Tdata(2:2:end,1))/2.^15;                      % 30720*12,PCI=0,root 25,subframe not align header,need sync
clear Tdata;


%% IQ data preprocess
% if (size(TimeData40sf,1) == 30720*40)
%     SF_num = 40; 
% elseif (size(TimeData40sf,1) == 30720*10)
%     SF_num = 10; 
% end;

if (size(TimeData40sf,1) == 30720*40) && (cs_para.TDdata_comb_flag == 1)
        TimeData12sf = TimeData40sf(30720*0+1:30720*0+30720*10,:) + TimeData40sf(30720*10+1:30720*10+30720*10,:) + TimeData40sf(30720*20+1:30720*20+30720*10,:) + TimeData40sf(30720*30+1:30720*30+30720*10,:);
        SF_num = 10;   
else
     TimeData12sf = TimeData40sf(30720*0+1:30720*0+30720*SF_num,:);
end;
clear TimeData40sf;



% TimeData12sf = timeDomainSig;

% figure(100);
% plot((real(TimeData12sf)));
% grid on;

figure(101);
plot((abs(TimeData12sf)));
grid on;


%% Downsampling 32 ratio
[TimeData,TimeDataDwsampling] = DwSampling_10_40ms(TimeData12sf,DelayT,SF_num);


% ut
% outdata = load('H:\ICS 数据 -4db 5000频偏\第二次频偏补偿后timing_signal_temp.txt');
% TimeDataDwsampling = outdata(:,1) + 1j*outdata(:,2);

%% SSS coarse timing search, the timing method1:SSS-PSSS
[Peak_sss_timing_value,  Peak_sss_timing_idx,noise_power_sss_coarse_timing] = Sss_coarse_timing_search(TimeDataDwsampling,DelayT);

     end;  % end of delayT,unit Test

% unit test
% Peak_sss_timing_idx = 4608;
     
%% first, FFO is estimated and compensation based CP_method
[TimeDataDwsampling_ffoCP, foe_cp_ffo] = Foe_cp_ffo_estimation(TimeData,TimeDataDwsampling,Peak_sss_timing_idx,SF_num);
if 0
    TimeDataDwsampling = TimeDataDwsampling_ffoCP;
end;


%% PSS detection

n_FFT = cs_para.N_FFT;
[du, fre_seq, time_seq] = PssGenAll(n_FFT);


 for iterN = 1:cs_para.iterNum



% then, joint IFO and pss peak search，now without IFO
for cfo_index = 0%-30:15:30
    cand_fo = cfo_index*1e3;
    comp_phase = exp(j*2*pi*cand_fo./(cs_para.N_FFT*sys_para.sub_carrier_spacing).*(0:size(TimeDataDwsampling,1)-1));
    TimeDataDwsamplingXcor = TimeDataDwsampling .*comp_phase.';
    
    
    peak = zeros(3,3);
    pos  = zeros(3,3);
    xcor  = zeros(3,30720/cs_para.Ratio*SF_num);
    Axcor  = zeros(3,30720/cs_para.Ratio*SF_num);

        
    for i = 1:3

        if (cs_para.split_xcorr_flag ==1)
            xcor(i,:) = circle_corr_split_square(TimeDataDwsamplingXcor.',time_seq(i,:),4);        % power calc
            Pxcor = xcor;
           
        else   % direct xcorr,no split segment
            xcor(i,:) = circle_corr(TimeDataDwsamplingXcor.',time_seq(i,:));
            Pxcor(i,:) = abs(xcor(i,:)).*abs(xcor(i,:));          
      
        end;   % xcorr power end
        
        if (cs_para.pss_slidewindow_xcorr_flag == 10)
            WLpss = 64;
            xcortmp(i,:) = [Pxcor(i,size(Pxcor,2)-WLpss/2+1:size(Pxcor,2)) Pxcor(i,:) Pxcor(i,1:WLpss/2)];
%             xcortmp(i,:) = [ones(1,WLpss/2) Pxcor(i,:) ones(1,WLpss/2)];
            for jpss = (WLpss/2)+1:size(xcor,2)+(WLpss/2)
                xcortmp2(i,jpss) = xcortmp(i,jpss)/mean(xcortmp(i,[jpss-WLpss/2:jpss-1 jpss+1:jpss+WLpss/2]),2);
%                 if (xcortmp2(i,jpss) == NaN)
%                     xcortmp2(i,jpss) = 0;
%                 end;
            end;
            
            xcor(i,:) = xcortmp2(i,[(WLpss/2)+1:size(xcor,2)+(WLpss/2)]);
            %                  figure();plot(xcor(i,:));grid on;
        elseif (cs_para.pss_slidewindow_xcorr_flag == 11)  
            WLpss = 64;
            xcortmp3 = zeros(3,960);
            xcortmp4 = [];
            for sf_index = 1:SF_num
                xcortmp(i,:) = [zeros(1,WLpss/2) Pxcor(i,30720/cs_para.Ratio*(sf_index-1)+1:30720/cs_para.Ratio*sf_index) zeros(1,WLpss/2)];
                for jpss = (WLpss/2)+1:size(xcortmp,2)-(WLpss/2)
                    xcortmp2(i,jpss) = xcortmp(i,jpss)/mean(xcortmp(i,[jpss-WLpss/2:jpss-1 jpss+1:jpss+WLpss/2]),2);
                    %                 if (xcortmp2(i,jpss) == NaN)
                    %                     xcortmp2(i,jpss) = 0;
                    %                 end;
                end;
                xcortmp3(i,:) = xcortmp2(i,[(WLpss/2)+1:size(xcortmp,2)-(WLpss/2)]);
                xcortmp4 = [xcortmp4 xcortmp3];
            end;
            
             xcor(i,:) = xcortmp4(i,:);
            
        end;  % silde window Pxcorr post process
        
        
        xcor_num22 = sum(reshape(xcor,3,4800,[]),3);             % comb 8 times for 40ms or 2times for 10ms;        
        Pxcor_num22(i,:) = abs(xcor_num22(i,:));      % search peak and calc noise
        P2xcor_num22 = Pxcor_num22;                   % plot and 
      
        
         % PSS threshold,   3 column for 3 Pss,  more line for more peaks;
         if (strcmp(cs_para.pss_noise_method,'PART'))
             [peak_pss(1,i), pos(1,i)] = max(Pxcor_num22(i,:));% 1st peak
             Pxcor_num22(i,pos(1,i)) = 0;
             [peak_pss(2,i), pos(2,i)] = max(Pxcor_num22(i,:));% 2st peak
             Pxcor_num22(i,pos(2,i)) = 0;
             noise_power_pss(i,1) = mean(Pxcor_num22(i,:));
             Ppss_th(i,1) = cs_para.pss_threshold*noise_power_pss(i,1);
         elseif(strcmp(cs_para.pss_noise_method,'ALL'))
             noise_power_pss(i,1) = mean(Pxcor_num22(i,:));
             Ppss_th(i,1) = cs_para.pss_threshold*noise_power_pss(i,1);
         end;
         
         %          markpeak = 1103;
         %          P2xcor_num22(:,[1:markpeak-8 markpeak+8:end]) = 0;
         
      
         figure((DelayT+cfo_index+100)*2+i);
         stem(P2xcor_num22(i,:));      grid on; hold on;
         plot([0 size(P2xcor_num22,2)],[Ppss_th(i,1) Ppss_th(i,1)],'r',[0 size(P2xcor_num22,2)],[noise_power_pss(i,1) noise_power_pss(i,1)],'c');
         ylabel('Pxcor combine');
         title(int2str(cand_fo));
         set (gcf,'Position',[i*380,100,400,300], 'color','w')
     
    end;   % end of 3 pss
    
%     
%     Peak_pss_timing_idx = 14;
%     pss_cand_WL = 10;
%     zoom_pss_xcor = P2xcor_num22(:,Peak_pss_timing_idx-pss_cand_WL/2:Peak_pss_timing_idx+pss_cand_WL/2);    
%     for i = 1:3
%         zoom_tmpnoise = sort(zoom_pss_xcor(i,:),'ascend');
%         zoom_noise_power_pss(i,1) = mean(zoom_tmpnoise(1,1:9));
%         zoom_Ppss_th(i,1) = cs_para.pss_threshold*zoom_noise_power_pss(i,1);
%         figure();
%          stem([1:size(zoom_pss_xcor,2)]+Peak_pss_timing_idx-6,zoom_pss_xcor(i,:),'b');      grid on; hold on;
%          plot([1+Peak_pss_timing_idx-6  size(zoom_pss_xcor,2)+Peak_pss_timing_idx-6],[zoom_Ppss_th(i,1) zoom_Ppss_th(i,1)],'r',[1+Peak_pss_timing_idx-6  size(zoom_pss_xcor,2)+Peak_pss_timing_idx-6],[zoom_noise_power_pss(i,1) zoom_noise_power_pss(i,1)],'c');
%          ylabel('zoom Pxcor combine');
%          title('zoom Pss Xcor power');
%          set (gcf,'Position',[i*380,100,400,300], 'color','w')
%     end;
%     
    
    
end;  % end of cand fo_index

% unit TEst
TimeDataDwsampling_reserve = TimeDataDwsampling;
TimeDataDwsampling = TimeDataDwsamplingXcor;


%% Pss reliability filter
%  cand_pss: 若干行*4列 ，4列分别是NID2No、pos、peakvalue、功率可靠峰flag，
%  后面在sss timing 后增加第5列定时可靠峰操作，此处暂不考虑不同ICFO下的pss峰值比较
%
[cand_pss_ca, cand_pss_ca_order] = Pss_reliability_filtering(P2xcor_num22,Ppss_th);
[cand_pss_ca_order2] = Pss_reliability_sss_timing_filtering(cand_pss_ca_order,TimeDataDwsampling,SF_num);  % 区分order与order2是为了便于回溯查看不满足sss timing的pss峰值点

    
    %% search Pss 3 peaks
    
    
    for i = 1:3
        [peak(i,1),pos(i,1)] = max(P2xcor_num22(i,:));                     % 4544/32+1=142+1=143; 30720/32 = 960;  960 + 143 =1103;
        %         [peak_ideal(i,1),pos_ideal(i,1)] = max(Axcor_ideal(i,:));   % true positon:160+2048+144+2048+144=4544,即Pss位置为4545； 30720+4545 = 35265;
    end;
    
    [peaktmp, NID2No] = max(peak(:,1));
     NID2 = NID2No-1;
   
    
   %%  coarse freq offset estimate and compensation, PSS,FCFO
   if (cs_para.unit_test_flag == 1)
      NID2No = 3;               %% Unit Test     NID2No = NID2 + 1;                %%%%%%%%%%%%%% Unit Test： 
      pos(NID2No,1) = 165; 
      NID2 = NID2No - 1;           % Unit Test
   end;
   
    if (cs_para.pss_freqoffset_comp_flag == 1)
        if (SF_num == 40)
            sumtmp = 0;
            for icom = 1:SF_num/5-1
                pos_i = pos(NID2No,1) + (icom -1)*30720/cs_para.Ratio*5;
                Tpss = TimeDataDwsampling(pos_i:pos_i + cs_para.N_FFT-1,1);
                r_pss = Tpss.*(time_seq(NID2No,:))';
                if strcmp(cs_para.pss_freqoffset_comp_method, 'SYS')
                    r_corse = sum(r_pss(33:64,1).*conj(r_pss(1:32,1)),1);
                    sumtmp = r_corse + sumtmp;
                elseif strcmp(cs_para.pss_freqoffset_comp_method, 'DSP')
                    rc = sum(r_pss(33:64,1)).*conj(sum(r_pss(1:32,1)));
                    sumtmp = rc + sumtmp;
                end;
            end;
            foe1 = angle(sumtmp)./pi*15000;
        else
            Tpss = TimeDataDwsampling(pos(NID2No,1):pos(NID2No,1)+ cs_para.N_FFT-1,1);
            r_pss = Tpss.*(time_seq(NID2No,:))';
            if strcmp(cs_para.pss_freqoffset_comp_method, 'SYS')
                r_corse = sum(r_pss(33:64,1).*conj(r_pss(1:32,1)),1);     %  bench method
                foe1 = angle(r_corse)./pi*15000;
            elseif strcmp(cs_para.pss_freqoffset_comp_method, 'DSP')
                rc = sum(r_pss(33:64,1)).*conj(sum(r_pss(1:32,1)));    % dsp method
                foe1 = angle(rc)./pi*15000;
            end;

        end;
        
%         close all;return;  % unit Test
        
        comp_phase = exp(-j*2*pi*foe1./(cs_para.N_FFT*sys_para.sub_carrier_spacing).*(0:size(TimeDataDwsampling,1)-1));
        TimeDataDwsampling2 = TimeDataDwsampling .*comp_phase.';
    else
        TimeDataDwsampling2 = TimeDataDwsampling;
    end;
   cs_para.pss_coarse_freq_offset = foe1;
   
    
    
    %% fine time offset pss cp, not valid
%     pos_up(NID2No,1) = (pos(NID2No,1)-1)*cs_para.Ratio + 1;
%     win_deltaT = [-cs_para.Ratio:cs_para.Ratio];
%     calc_len = 80;
%     for n_Fine = 1:size(win_deltaT,2)
%         xcor_cp(n_Fine,1) = (TimeData(pos_up(NID2No,1)+win_deltaT(n_Fine)-calc_len:pos_up(NID2No,1)+win_deltaT(n_Fine)-1,1)).'*...
%                             conj(TimeData(pos_up(NID2No,1)+win_deltaT(n_Fine)-calc_len+sys_para.N_FFT:pos_up(NID2No,1)+win_deltaT(n_Fine)-1+sys_para.N_FFT,1));
%         Pxcor_cp(n_Fine,1) = abs(xcor_cp(n_Fine,1));                
%     end;
%     figure();plot(Pxcor_cp,'-o');grid on;
%     
    
   
    %% % extract SSS, 30720-2048+1 = 28673;  960-64+1=897; sf0 sss
    if (cs_para.sss_detec_halfframe == 0)   
        Pos_pss1 = pos(NID2No,1);
    else
        Pos_pss1 = pos(NID2No,1)+4800;
    end;
    
%     NID2 = NID2No - 1;           % Unit Test
    
% for NID2 = 0:2
    window = [0];     % TA offset;
    if (cs_para.Mode == 'TDD')
        Pos_sss_origin = Pos_pss1 - 4544/32 - 2048/32;
    else 
        Pos_sss_origin = Pos_pss1 - ceil(144/32) - 2048/32;
    end;
    
    if (Pos_sss_origin < 0)
        Pos_sss_origin = Pos_sss_origin + 30720*5/32 ;
    end;
    
    [Lfreq_sss_all, Ltime_sss_all] = SssGenAll(n_FFT,NID2);    

    
    for n_FineTiming = 1:size(window,2)
        %         Pos_sss1 = Pos_sss_origin + window(n_FineTiming);
        Pxcorr_sss_f_tmp = zeros(336,1);
        Pxcorr_sss_f_tmp2 = zeros(336,1);
        if (cs_para.sss_xcorr_comb_flag == 1)%&&(SF_num == 40)
            for iss = 1:SF_num/10
                
                Pos_sss1 = Pos_sss_origin + window(n_FineTiming) + (iss-1)*30720*10/cs_para.Ratio;
                
                %                 for iss2 = 1:2
                %                     if (iss2 == 2)
                %                         Pos_sss1 = Pos_sss1 + 4800;
                %                     end;
                
                
                Rt_sss = TimeDataDwsampling2(Pos_sss1:Pos_sss1+63,1);
                Rf_sss = sqrt(1.0/n_FFT)*fft(Rt_sss,64,1);
                zc_sss = [Rf_sss(34:64,1);Rf_sss(2:32,1)];
                %
                
                if (cs_para.sss_equal_flag == 1)
                    Pos_pss_tmp = Pos_pss1 + (iss-1)*30720*10/cs_para.Ratio;
                    Rt_pss = TimeDataDwsampling2(Pos_pss_tmp:Pos_pss_tmp+63,1);
                    Rf_pss = sqrt(1.0/n_FFT)*fft(Rt_pss,64,1);
                    Rzc_pss = [Rf_pss(34:64,1);Rf_pss(2:32,1)];
                    
                    Hls_pss = Rzc_pss.'./du(NID2+1,:);
                    %                                figure();plot(real(Hls_pss));grid on;
                    
                    if strcmp(cs_para.sss_equal_pss_Hpost_method,'TD')
                        % alt1:Hls reduce noise; tap = 0;
                        hls_pss = sqrt(n_FFT)*ifft(Hls_pss,62);
                        %                                 figure();plot(abs(hls_pss));grid on;
                        hls_pss(1,4:62) = 0;
                        H_pss = sqrt(1.0/n_FFT)*fft(hls_pss);
                        %                            figure();plot(real(H_pss));grid on;
                        
                        
                    elseif strcmp(cs_para.sss_equal_pss_Hpost_method,'FILTER')
                        % alt2:filter pss sss hls
                        WL = 15;
                        for i = 1:floor(WL/2)
                            Hpss(1,i) = sum(Hls_pss(1,1:WL),2)/WL;
                        end;
                        for i = floor(WL/2)+1:size(Hls_pss,2)-ceil(WL/2)+1
                            Hpss(1,i) = sum(Hls_pss(1,i-floor(WL/2):i+ceil(WL/2)-1),2)/WL;
                        end;
                        for i = size(Hls_pss,2)-ceil(WL/2)+2:size(Hls_pss,2)
                            Hpss(1,i) = sum(Hls_pss(1,size(Hls_pss,2)-WL+1:size(Hls_pss,2)),2)/WL;
                        end;
                        
                        H_pss = Hpss;
                        %                          figure();plot(abs(Hls_pss));grid on;
                        %                                    figure();plot(real(H_pss));grid on;
                        
                    end;
                    
                    %                         eq_sss = inv(eye(62).*(H_pss'*H_pss))*H_pss'.*zc_sss;     % use pss che to equal sss
                    eq_sss = H_pss'.*zc_sss;
                    %                         eq_sss = H_pss'.*zc_sss/norm(H_pss)*sqrt(62);
                    zc_sss = eq_sss;
                    
                end;   % finish eq sss
                
                eq_zc_sss(:,iss) = zc_sss;
            end;        % end of comb 40 subframe
            
            % xcorr_sss_t = Ltime_sss_all*conj(Rt_sss);
            % Axcorr_sss_t = abs(xcorr_sss_t);
            %
            % figure();
            % stem(Axcorr_sss_t(:,1));      grid on;
            % ylabel('Axcor_SSS_time');
            
            if (cs_para.sss_equal_data_realize_flag == 1)
                eq_zc_sss = real(eq_zc_sss);
            end;
            
            if strcmp(cs_para.sss_xcorr_calc_method,'336full')
                xcorr_sss_f = Lfreq_sss_all*conj(eq_zc_sss);            % 336*62   *   62*4/1
            else strcmp(cs_para.sss_xcorr_calc_method,'m0m1')
                [c0 c1] = c0c1Gen(NID2);
                zc_sss_even = eq_zc_sss(1:2:end,:).*repmat(c0.',1,size(eq_zc_sss,2));                              
                [sbase] = sbaseGen;
                for comb_index = 1:size(zc_sss_even,2)                    
                    SSS_xorr_even(:,comb_index) = circle_corr(sbase,zc_sss_even(:,comb_index).');
                end;                
                xcorr_sss_f = SSS_xorr_even;                
            end;
            
            %                     if (iss2 == 2)
            %                         xcorr_sss_f = [Lfreq_sss_all(169:336,:);Lfreq_sss_all(1:168,:)]*conj(zc_sss);
            %                     end;
            
            if strcmp(cs_para.sss_xcorr_comb_method,'VECTOR')
                Pxcorr_sss_f_tmp2 = sum(xcorr_sss_f,2);% + Pxcorr_sss_f_tmp2;  % xcorr comb
                Pxcorr_sss_f_tmp = abs(Pxcorr_sss_f_tmp2).*abs(Pxcorr_sss_f_tmp2);
            else strcmp(cs_para.sss_xcorr_comb_method,'POWER')
                Pxcorr_sss_f_tmp = sum((abs(xcorr_sss_f).*abs(xcorr_sss_f)),2);% + Pxcorr_sss_f_tmp;  % power comb
            end;
                    
%                 end;      % end of comb sf0 and sf5;
           
            
%         else        % only 10ms, only sf0 or sf5,no comb
%             Pos_sss1 = Pos_sss_origin + window(n_FineTiming);
%             Rt_sss = TimeDataDwsampling2(Pos_sss1:Pos_sss1+63,1);
%             Rf_sss = sqrt(1.0/n_FFT)*fft(Rt_sss,64,1);
%             zc_sss = [Rf_sss(34:64,1);Rf_sss(2:32,1)];
% %             zc_sss(1:2:end,1) = zc_sss(1:2:end,1).*c0.';
% %             zc_sss(2:2:end,1) = zc_sss(2:2:end,1).*c1.';
%                     
%             if (cs_para.sss_equal_flag == 1)
%                 Rt_pss = TimeDataDwsampling2(Pos_pss1:Pos_pss1+63,1);
%                 Rf_pss = sqrt(1.0/n_FFT)*fft(Rt_pss,64,1);
%                 Rzc_pss = [Rf_pss(34:64,1);Rf_pss(2:32,1)];
%                 
%                 Hls_pss = Rzc_pss.'./du(NID2+1,:);
%                 hls_pss = sqrt(n_FFT)*ifft(Hls_pss);
%                 %                     figure();plot(real(Hls_pss));grid on;
%                 %                     figure();plot(abs(hls_pss));grid on;
%                 
%                 
%                 
%                 if strcmp(cs_para.sss_equal_pss_Hpost_method,'TD')
%                     % alt1:Hls reduce noise; tap = 0;
%                     hls_pss(1,4:60) = 0;
%                     H_pss = sqrt(1.0/n_FFT)*fft(hls_pss);
%                     %                         figure();plot(real(H_pss));grid on;
%                     
%                     
%                 else
%                     % alt2:filter pss sss hls
%                     WL = 15;
%                     for i = 1:floor(WL/2)
%                         Hpss(1,i) = sum(Hls_pss(1,1:WL),2)/WL;
%                     end;
%                     for i = floor(WL/2)+1:size(Hls_pss,2)-ceil(WL/2)+1
%                         Hpss(1,i) = sum(Hls_pss(1,i-floor(WL/2):i+ceil(WL/2)-1),2)/WL;
%                     end;
%                     for i = size(Hls_pss,2)-ceil(WL/2)+2:size(Hls_pss,2)
%                         Hpss(1,i) = sum(Hls_pss(1,size(Hls_pss,2)-WL+1:size(Hls_pss,2)),2)/WL;
%                     end;
%                     
%                     H_pss = Hpss;
%                 end;
%                 
% %                 eq_sss = inv(eye(62).*(H_pss'*H_pss))*H_pss'.*zc_sss;       % use pss che to equal sss
%                 eq_sss = H_pss'.*zc_sss;
%                 zc_sss = eq_sss;
%                 
%             end;       % end sss eqal calc
%             
%             xcorr_sss_f = Lfreq_sss_all*conj(zc_sss);
%             Pxcorr_sss_f_tmp = abs(xcorr_sss_f).*abs(xcorr_sss_f);  % power;
        end;      % end 40ms comb xcorr calc
        
        Pxcorr_sss_f(:,n_FineTiming) = Pxcorr_sss_f_tmp;
        P2xcorr_sss_f = Pxcorr_sss_f;      % for plot and search
        
        if strcmp(cs_para.sss_noise_method, 'PART' )
            [peak_sss(1,n_FineTiming), tmpNID1No(1,n_FineTiming)] = max(Pxcorr_sss_f(:,n_FineTiming));% 1st peak
            Pxcorr_sss_f(tmpNID1No(1,n_FineTiming),n_FineTiming) = 0;
            [peak_sss(2,n_FineTiming), tmpNID1No(2,n_FineTiming)] = max(Pxcorr_sss_f(:,n_FineTiming)); % 2nd peak
            Pxcorr_sss_f(tmpNID1No(1,n_FineTiming),n_FineTiming) = 0;
            noise_power(1,n_FineTiming) = mean(Pxcorr_sss_f(:,n_FineTiming));
            Psss_th(1,n_FineTiming) = cs_para.sss_threshold*noise_power(1,n_FineTiming);
        elseif strcmp(cs_para.sss_noise_method,'ALL' )
            noise_power(1,n_FineTiming) = mean(P2xcorr_sss_f(:,n_FineTiming));
            Psss_th(1,n_FineTiming) = cs_para.sss_threshold*noise_power(1,n_FineTiming);
            
        end;           % thresold calc,only for m0
        
        figure();
        plot(P2xcorr_sss_f(:,n_FineTiming));      grid on;
        hold on; plot([0 size(Pxcorr_sss_f,1)],[Psss_th(1,n_FineTiming) Psss_th(1,n_FineTiming)],'r',[0 size(Pxcorr_sss_f,1)],[noise_power(1,n_FineTiming) noise_power(1,n_FineTiming)],'c');
        ylabel('Pxcor_SSS_freq');
        title('sss or 1st mi');
        set (gcf,'Position',[380,500,400,300], 'color','w');
        
        if strcmp(cs_para.sss_xcorr_calc_method,'336full')
            % search peak SSS and NID1
            [peakValueSss, NoInWin] = max(peak_sss(1,:));       % search which position in window
            [peakValueSss, NID1No] = max(P2xcorr_sss_f(:,NoInWin));% search NID1
            NID1 = NID1No-1;
            
            if (NID1No > 168)
                NID1 = NID1No-1-168;
            end;
            
            % pci
            cs_para.detPCI = [  NID1*3+NID2];
            
        else strcmp(cs_para.sss_xcorr_calc_method,'m0m1')
            mi_even_cad = tmpNID1No(:,1);
            for sss_dec_index = 1:cs_para.sss_xcorr_peak_select_num                  
                m0index = mi_even_cad(sss_dec_index,1)-1;
                [z1mi] = z1miGen(m0index);
                zc_sss_odd = eq_zc_sss(2:2:end,:).*repmat(c1.',1,size(eq_zc_sss,2)).*repmat(z1mi.',1,size(eq_zc_sss,2));
                for comb_index = 1:size(zc_sss_odd,2)                    
                    SSS_xorr_odd(:,comb_index) = circle_corr(sbase,zc_sss_odd(:,comb_index).');
                end;
                
                if strcmp(cs_para.sss_xcorr_comb_method,'VECTOR')
                    Pxcorr_sss_odd_tmp2 = sum(SSS_xorr_odd,2);% + Pxcorr_sss_f_tmp2;  % xcorr comb
                    Pxcorr_sss_odd_tmp = abs(Pxcorr_sss_odd_tmp2).*abs(Pxcorr_sss_odd_tmp2);
                else strcmp(cs_para.sss_xcorr_comb_method,'POWER')
                    Pxcorr_sss_odd_tmp = sum((abs(SSS_xorr_odd).*abs(SSS_xorr_odd)),2);% + Pxcorr_sss_f_tmp;  % power comb
                end;
                
                figure();
                plot(Pxcorr_sss_odd_tmp);      grid on;
                hold on; plot([0 size(Pxcorr_sss_odd_tmp,1)],[Psss_th(1,n_FineTiming) Psss_th(1,n_FineTiming)],'r',[0 size(Pxcorr_sss_odd_tmp,1)],[noise_power(1,n_FineTiming) noise_power(1,n_FineTiming)],'c');
                ylabel('Pxcor_SSS_m1');
                title('2nd mi'); set (gcf,'Position',[760,500,400,300], 'color','w');
                
                [m1value,m1index] = max(Pxcorr_sss_odd_tmp);
                m1index = m1index - 1;
                
                if (m0index < m1index)
                    m0 = m0index;
                    m1 = m1index;
                    iHalf = 0;
                else
                    m0 = m1index;
                    m1 = m0index;
                    iHalf = 1;
                end;
                
                bSuccess = true;
                if ((m1-m0>7) || (m1 == m0) ||(m1-m0 == 7 && m0>2))
                    bSuccess = false;                    
                end;
                if (bSuccess)
                    m_tmp = 31*(m1-m0-1)+m0;
                    k_tmp = floor(m_tmp/30);
                    NID1 = m_tmp - floor(k_tmp*(k_tmp+1)/2);
                    tmpPCI = 3*NID1+NID2;
                else
                    tmpPCI = -1;
%                     return;
                end;
                if (iHalf == 0)
                    NID1No = NID1 + 1;
                else
                    NID1No = NID1 + 1 +168;
                end;
                
                cs_para.detPCI = [cs_para.detPCI tmpPCI];
            end;   % end of pci dec num
        end;     % end of sss dec method
        
    end;    % end of n_FineTiming
 
%      end;    % end of freqoffset,unit test;
   
% end;    % end of NID2 loop
    
%     if(peakValueSss < Psss_th(1,NID1No))
%         NID1 = -1;                      % invalid NID1;
%         continue;
%     end;
 
%% ICFO,PSS， use origin data with pss coarse freq offset compensation
Pos_pss_origin_up = (Pos_pss1-1)*cs_para.Ratio+1;
Pxcor_fd_pss_sum = 0;
for iss = 1:SF_num/5
    Pos_pss_up = Pos_pss_origin_up + (iss-1)*153600;
    Rpss_up = TimeData(Pos_pss_up:Pos_pss_up+2048-1,1);
    nn = 1:2048;
    foe_pss_coarse_com = exp(-j*2*pi*foe1./(2048*15000).*(nn-1));
     Rpss_up = Rpss_up.*foe_pss_coarse_com.';
    phy_sc_index = phy_scmap(sys_para.N_FFT,sys_para.N_RB_DL,sys_para.Nsc_RB,'OFDMA');
    Rpss_sf = ofdm_fft(Rpss_up,sys_para.N_FFT,1);
    Rpss_sf_FD   = Rpss_sf([phy_sc_index],:);
    xcorr_pss_FD = circle_corr(Rpss_sf_FD.',du(NID2No,:));
    Pxcor_fd_pss = abs(xcorr_pss_FD).*abs(xcorr_pss_FD);
    Pxcor_fd_pss_sum = Pxcor_fd_pss_sum + Pxcor_fd_pss;
end;
noise_power_FD_pss = mean(Pxcor_fd_pss_sum);
Ppss_ICO_th = cs_para.pss_ICO_threshold*noise_power_FD_pss;
figure();
plot([1:1200],Pxcor_fd_pss_sum,'b',[0 1200],[Ppss_ICO_th Ppss_ICO_th],'r',[0 1200],[noise_power_FD_pss noise_power_FD_pss],'c');grid on;
ylabel('Pxcor_PSS_ICO');

%% ICFO,SSS
Pos_sss_origin_up = (Pos_sss_origin-1)*cs_para.Ratio+1;
 Pxcor_fd_sss_sum = 0;
for iss = 1:SF_num/10
    Pos_sss_up = Pos_sss_origin_up + (iss-1)*307200;
    Rsss_up = TimeData(Pos_sss_up:Pos_sss_up+2048-1,1);
    nn = 1:2048;
    foe_sss_coarse_com = exp(-j*2*pi*foe1./(2048*15000).*(nn-1));
    Rsss_up = Rsss_up.*foe_sss_coarse_com.';
    phy_sc_index = phy_scmap(sys_para.N_FFT,sys_para.N_RB_DL,sys_para.Nsc_RB,'OFDMA');
    Rsss_sf = ofdm_fft(Rsss_up,sys_para.N_FFT,1);
    Rsss_sf_FD   = Rsss_sf([phy_sc_index],:);
    xcorr_sss_FD = circle_corr(Rsss_sf_FD.',Lfreq_sss_all(NID1No,:));
    Pxcor_fd_sss = abs(xcorr_sss_FD).*abs(xcorr_sss_FD);
    Pxcor_fd_sss_sum = Pxcor_fd_sss_sum + Pxcor_fd_sss;
end;
noise_power_FD_sss = mean(Pxcor_fd_sss_sum);
Psss_ICO_th = cs_para.sss_ICO_threshold*noise_power_FD_sss;
figure();
plot([1:1200],Pxcor_fd_sss_sum,'b',[0 1200],[Psss_ICO_th Psss_ICO_th],'r',[0 1200],[noise_power_FD_sss noise_power_FD_sss],'c');grid on;
ylabel('Pxcor_SSS_ICO');

    
%% %  RSRP RSRQ RSSI , sf0 calc, for ver dsp
%      NID2No = 2;                                 % UT,  NID2No = NID2 + 1;   1 2 3 for 0 1 2
%      pos(NID2No,1) = 1104;                       % UT
if (cs_para.unit_test_flag == 1)
    cs_para.detPCI = 353;    
end; % UT


if (cs_para.M1_RPRQ_flag == 1)
    
    if (NID1No > 168)                        % peak is sf 6, the next is sf 0;
        if (cs_para.Mode == 'TDD')
            Msf_start =  (pos(NID2No,1)-1)*32+1-144*2-160-2048*2+4*30720;
        else
            Msf_start =  (pos(NID2No,1)-1)*32+1-144*6-160-2048*6+5*30720;
        end;
    else
        if (cs_para.Mode == 'TDD')
            if (pos(NID2No,1) > 1102)             % peak is sf 1, and the before sf is full sf0;
                Msf_start =  (pos(NID2No,1)-1)*32+1-144*2-160-2048*2-30720;
            else                                   % peak is sf 1, and the before sf is part sf0;
                Msf_start =  (pos(NID2No,1)-1)*32+1-144*2-160-2048*2+9*30720;
            end;
        else
            if (pos(NID2No,1) > 416)             % peak is sf 0, and is full sf0;
                Msf_start =  (pos(NID2No,1)-1)*32+1-144*6-160-2048*6;
            else                                   % peak is sf 0, and is part sf0;
                Msf_start =  (pos(NID2No,1)-1)*32+1-144*6-160-2048*6+10*30720;
            end;;
        end;
    end;
    
    %        cs_para.detPCI = 398;    %;                % unit test, assign PCI;
%            Msf_start = Msf_start + 30720*7 - 624;                                % unit test, time offset
    %
    TimeDataDouble = [TimeData;TimeData];           % for extract sss,need enough data size;
    Msf = TimeDataDouble(Msf_start:Msf_start+30720-1,1); clear TimeDataDouble;
    % Msf = [TimeData(Msf_start:Msf_start+30720-1-624,1);zeros(624,1)];
        
    nn = 1:30720;
    foe_sf0_coarse_com = exp(-j*2*pi*foe1./(2048*15000).*(nn-1));
    Msf = Msf.*foe_sf0_coarse_com.';  
    
    Msf_t = Msf.';                                     % Rx*sampling_num;
    [M_TD M_TD0 M_FD M_FD0 UT_FD_a] = DL_ofdm_demodulation(Msf_t,Msf_t);
    
    figure();
    plot(abs(reshape(UT_FD_a,[],1)));
    grid on;
    
    
    sf0_td = squeeze(M_TD);
    sf0_fd = squeeze(M_FD);
    figure();
    scatter(real(sf0_fd([570:631],6)),imag(sf0_fd([570:631],6)));   % SSS Constellation，
    grid on;    title('SSS Constellation');
    
    figure();
    scatter(real(sf0_fd([570:631],7)),imag(sf0_fd([570:631],7)));   % SSS Constellation，
    grid on;    title('PSS Constellation');
    
    % write sf0 TD FD data;
    if(TestVectorPrintFlag)
        fid3 = fopen('ICS_SF2_TD.dat','w+');
        for i = 1:1:size(Msf)
            fprintf(fid3,'%10.6f     %10.6f \n' ,real(Msf(i,1)),imag(Msf(i,1)));
        end;
        
        CRS_FD_Print = reshape(squeeze(M_FD),[],1);
        fid4 = fopen('ICS_SF0_FD.dat','w+');
        for i = 1:1:size(CRS_FD_Print)
            fprintf(fid3,'%10.6f     %10.6f \n' ,real(CRS_FD_Print(i,1)),imag(CRS_FD_Print(i,1)));
        end;
        fclose('all');
    end;
    
    %     UT_FD = squeeze(UT_FD_a(1,:,:));
    %     UT_FD_plot = reshape(UT_FD_a,[],1);
    %
    %     figure();
    %     semilogy((abs(UT_FD_plot)));
    %     grid on;
    %     figure();
    %     plot((imag(UT_FD_plot)));
    %     grid on;
    
    CRS_symbol = [0 4 7 11];
    CRS_TD = squeeze(M_TD(1,:,CRS_symbol+1));
    CRS_FD = squeeze(M_FD(1,:,CRS_symbol+1));
        
    % RSSI, CRS 0 4 7 11
    Pcrs_td = abs(CRS_TD).*abs(CRS_TD);                              % method 1: TD
    RSSI = 1/size(CRS_TD,2)*1/size(CRS_TD,1)*sum(sum(Pcrs_td,1),2);
    Prssi_fd = abs(CRS_FD).*abs(CRS_FD);                             % method 2: FD
    RSSI_FD = 1200/2048*1/size(CRS_FD,2)*1/size(CRS_FD,1)*sum(sum(Prssi_fd,1),2);
    RSSI_dBv_fxd = 10*log10(RSSI);
    RSSI_dBv_flp = RSSI_dBv_fxd;% - 10*log10(2.^22);
    RSSI_dBm = RSSI_dBv_flp - 4.1 + 13 - LNAGain - AD9361Gain;
    
    
    % RSSI，pss，sss
    pss_start_up =  Msf_start + 144*2+160+2048*2+30720;%(pos(NID2No,1)-1)*32+1;
    sss_start_up =  pss_start_up-144*2-160-2048*3;
    pss_up = TimeData(pss_start_up:pss_start_up+2048-1,1); 
    sss_up = TimeData(sss_start_up:sss_start_up+2048-1,1); 
    Ppss_td = abs(pss_up).*abs(pss_up);   
    Psss_td = abs(sss_up).*abs(sss_up);   
    PSS_SI = 1/size(Ppss_td,1)*sum(Ppss_td,1);
    SSS_SI = 1/size(Psss_td,1)*sum(Psss_td,1);
    PSS_SI_dBv_flp = 10*log10(PSS_SI);
    SSS_SI_dBv_flp = 10*log10(SSS_SI);
    PSS_SI_dBm = PSS_SI_dBv_flp - 4.1 + 13 - LNAGain - AD9361Gain;
    SSS_SI_dBm = SSS_SI_dBv_flp - 4.1 + 13 - LNAGain - AD9361Gain;
    cs_para.det_PSS_SI = [PSS_SI_dBv_flp PSS_SI_dBm];
    cs_para.det_SSS_SI = [SSS_SI_dBv_flp SSS_SI_dBm];  
    
    
    
    % RSRP
    
    vshift = mod(cs_para.detPCI,6);
    vshift1 = mod(vshift + 3,6);
    CRS07 = CRS_FD(vshift+1:6:end,[1 3 ]);
    CRS411 = CRS_FD(vshift1+1:6:end,[2 4 ]);
    Rcrs = [CRS07(:,1) CRS411(:,1) CRS07(:,2) CRS411(:,2)];
    Pcrs_fd = abs(Rcrs).*abs(Rcrs);
    
    figure();
    scatter(real(reshape(Rcrs,[],1)),imag(reshape(Rcrs,[],1)));   % CRS Constellation，
    grid on;    title('CRS Constellation');
    set (gcf,'Position',[400,100,400,300], 'color','w')
    
    Pt = 1/2048*1/size(Pcrs_fd,1)*1/size(Pcrs_fd,2)*sum(sum(Pcrs_fd,1),2);
    if strcmp(cs_para.M1_method,'DSP')
        RSRP = Pt;
    elseif strcmp(cs_para.M1_method,'SYS')
        CRS_table = load('crs_matrix_normalCP.mat');
        CRS_local = squeeze(CRS_table.crs_matrix(cs_para.detPCI(1)+1,:,[1 2 3 4])); clear CRS_table;
        H_CRS = Rcrs./CRS_local;
        sub_H = [H_CRS(:,3)-H_CRS(:,1);H_CRS(:,4)-H_CRS(:,2)];
        Pn = abs(sub_H).*abs(sub_H);
        Pn_ave = 1/2048*1/2*1/size(Pn,1)*sum(Pn,1);
        RSRP = Pt - Pn_ave;
    end;
    RSRP_dBv_fxd = 10*log10(RSRP);
    RSRP_dBv_flp = RSRP_dBv_fxd;% -10*log10(2.^22);
    RSRP_dBm = RSRP_dBv_flp  - 4.1 + 13 - LNAGain - AD9361Gain;
    
    % RSRQ
    RSRQ_dB = RSRP_dBm - RSSI_dBm + 10*log10(sys_para.N_RB_DL);
    
    cs_para.det_RSRP = [RSRP_dBv_flp RSRP_dBm];
    cs_para.det_RSSI = [RSSI_dBv_flp RSSI_dBm];
    cs_para.det_RSRQ = RSRQ_dB;
    
 end;   % end of M1

%% CRS CTO and CFO estimation and compensation； FCFO within -1k and 1k
[foe_crs to_crs] = TimingFreqOffset(H_CRS);
cs_para.crs_freq_offset = foe_crs;
cs_para.crs_timing_offset = to_crs;


    
    %% Pss sf0/5 and SSS sf0 interfere cancel;
%      NID2No = 2; %UT
%      NID1Notmp = [119 133 133];
     NID1Notmp(1,NID2No) = NID1No;
%      pos(NID2No,1) = 3888;  %UT
%      for NID2 = 0:2
         NID2No = NID2 + 1;
         for iss = 1:SF_num/10
             Pos_sss1 = Pos_sss_origin + (iss-1)*30720*10/cs_para.Ratio;
             I_tmp = zeros(size(TimeDataDwsampling,1),1);
             for In = 1:2
                 if (In ==1)
                     Pos_sss1 = Pos_sss1;
                     NID1Notmp2 = NID1Notmp(1,NID2No);
                 else
                     Pos_sss1 = Pos_sss1 + 4800;
                     if (NID1Notmp(1,NID2No) > 168)
                         NID1Notmp2 = NID1Notmp(1,NID2No) - 168;
                     elseif (NID1Notmp(1,NID2No) < 169)
                         NID1Notmp2 = NID1Notmp(1,NID2No) + 168;
                     end;
                 end;
                 
                 Pos_pss1 = Pos_sss1 + 4544/cs_para.Ratio + 2048/cs_para.Ratio;
                 
                 Rt_pss = TimeDataDwsampling2(Pos_pss1:Pos_pss1+63,1);
                 Rf_pss = sqrt(1.0/n_FFT)*fft(Rt_pss,64,1);
                 zc_pss = [Rf_pss(34:64,1);Rf_pss(2:32,1)];
                 
                 Hls_pss = zc_pss.'./du(NID2No,:);
                 hls_pss = sqrt(n_FFT)*ifft(Hls_pss);
%                  figure();plot(real(Hls_pss));grid on;
%                  figure();plot(abs(hls_pss));grid on;
                 
                 Rt_sss = TimeDataDwsampling2(Pos_sss1:Pos_sss1+63,1);
                 Rf_sss = sqrt(1.0/n_FFT)*fft(Rt_sss,64,1);
                 zc_sss = [Rf_sss(34:64,1);Rf_sss(2:32,1)];
                 Hls_sss = zc_sss.'./Lfreq_sss_all(NID1Notmp2,:);
                 hls_sss = sqrt(n_FFT)*ifft(Hls_sss);
%                  figure();plot(real(Hls_sss));grid on;
%                  figure();plot(abs(hls_sss));grid on;
                 
                 if strcmp(cs_para.H_post_method,'TD')
                     % alt1:Hls reduce noise; tap = 0;
                     hls_pss(1,4:60) = 0;
                     H_pss = sqrt(1.0/n_FFT)*fft(hls_pss);
%                      figure();plot(real(H_pss));grid on;
                     
                     hls_sss(1,4:60) = 0;
                     H_sss = sqrt(1.0/n_FFT)*fft(hls_sss);
%                      figure();plot(real(H_sss));grid on;
                 elseif strcmp(cs_para.H_post_method,'FD')
                     % alt2:filter pss sss hls
                     hcoee = [  -0.07330612120138,  0.01719914911156,  0.08241779197212,   0.1681761005048,...
                         0.2408883292293,   0.2692821843831,   0.2408883292293,   0.1681761005048,...
                         0.08241779197212,  0.01719914911156, -0.07330612120138];
                     Htmp = conv(Hls_pss,hcoee);
                     H_pss = Htmp(6:67);
                     Htmp = conv(Hls_sss,hcoee);
                     H_sss = Htmp(6:67);
                 elseif strcmp(cs_para.H_post_method,'FILTER')
                     WL = 15;
                     for i = 1:floor(WL/2)
                         Hpss(1,i) = sum(Hls_pss(1,1:WL),2)/WL;
                     end;
                     for i = floor(WL/2)+1:size(Hls_pss,2)-ceil(WL/2)+1
                         Hpss(1,i) = sum(Hls_pss(1,i-floor(WL/2):i+ceil(WL/2)-1),2)/WL;
                     end;
                     for i = size(Hls_pss,2)-ceil(WL/2)+2:size(Hls_pss,2)
                         Hpss(1,i) = sum(Hls_pss(1,size(Hls_pss,2)-WL+1:size(Hls_pss,2)),2)/WL;
                     end;
                     
                     for i = 1:floor(WL/2)
                         Hsss(1,i) = sum(Hls_sss(1,1:WL),2)/WL;
                     end;
                     for i = floor(WL/2)+1:size(Hls_sss,2)-ceil(WL/2)+1
                         Hsss(1,i) = sum(Hls_sss(1,i-floor(WL/2):i+ceil(WL/2)-1),2)/WL;
                     end;
                     for i = size(Hls_sss,2)-ceil(WL/2)+2:size(Hls_sss,2)
                         Hsss(1,i) = sum(Hls_sss(1,size(Hls_sss,2)-WL+1:size(Hls_sss,2)),2)/WL;
                     end;
                     
                     H_pss = Hpss;
                     H_sss = Hsss;
                     
                     Pss_RP = 10*log10(sum(abs(H_pss).*abs(H_pss))/size(H_pss,2)/2048);
                     Sss_RP = 10*log10(sum(abs(H_sss).*abs(H_sss))/size(H_sss,2)/2048);
                     PSS_RP_dBm = Pss_RP  - 4.1 + 13 - LNAGain - AD9361Gain;
                     SSS_RP_dBm = Sss_RP  - 4.1 + 13 - LNAGain - AD9361Gain;
                     cs_para.det_PssRP = [Pss_RP PSS_RP_dBm];
                     cs_para.det_SssRP = [Sss_RP SSS_RP_dBm];
                     Pss_RQ_dB = PSS_RP_dBm - PSS_SI_dBm + 10*log10(sys_para.N_RB_DL);
                     Sss_RQ_dB = SSS_RP_dBm - PSS_SI_dBm + 10*log10(sys_para.N_RB_DL);
                     cs_para.det_Pss_RQ = Pss_RQ_dB;
                     cs_para.det_Sss_RQ = Sss_RQ_dB;
                     
                 end;
                 
                 % restruct interfere PSS SSS
                 if strcmp(cs_para.inf_restr_che_source, 'PSS')
                     Hif = H_pss;
                 elseif strcmp(cs_para.inf_restr_che_source, 'SSS')
                     Hif = H_sss;
                 end;
                 If_pss1 = du(NID2No,:).*Hif;     % use sss che to restruct pss
                 If_pss1 = If_pss1.';
                 I_seq = zeros(n_FFT,1);
                 I_seq(2:32) = If_pss1(32:62);    %  第1个为DC
                 I_seq(-30+n_FFT:n_FFT) = If_pss1(1:31);
                 It = sqrt(n_FFT)*ifft(I_seq);
                 
                 I_tmp(Pos_pss1:Pos_pss1+63,1) =  It;
                 
                 if strcmp(cs_para.inf_restr_che_source, 'PSS')
                     Hif = H_pss;
                 elseif strcmp(cs_para.inf_restr_che_source, 'SSS')
                     Hif = H_sss;
                 end;
                 If_sss1 = Lfreq_sss_all(NID1Notmp2,:).*Hif;
                 If_sss1 = If_sss1.';
                 I_seq = zeros(n_FFT,1);
                 I_seq(2:32) = If_sss1(32:62);    %  第1个为DC
                 I_seq(-30+n_FFT:n_FFT) = If_sss1(1:31);
                 It = sqrt(n_FFT)*ifft(I_seq);
                 I_tmp(Pos_sss1:Pos_sss1+63,1) =  It;
                 
                 if (cs_para.origin_data_reserve_flag == 0 )
                     IF_comp_phase = exp(j*2*pi*foe1./(cs_para.N_FFT*sys_para.sub_carrier_spacing).*(0:size(TimeDataDwsampling,1)-1));
                     I_tmp = I_tmp .*IF_comp_phase.';
                 end;
                 
             end;  % end In=1 2. that is 10ms
             
             
             % substrct interfere PSS and SSS
             if (cs_para.origin_data_reserve_flag == 0 )
                 Time2 = TimeDataDwsampling - I_tmp;
                 TimeDataDwsampling = Time2;                 
             elseif (cs_para.origin_data_reserve_flag == 1)
                 Time2 = TimeDataDwsampling2 - I_tmp;
                 TimeDataDwsampling = Time2;
             end;
                          
         end;    % end of 40ms
%      end;   % end of NID2

%          pause();

close all;
 end;            % end iterative;
