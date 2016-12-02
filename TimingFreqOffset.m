% H_ls:  200*4 in sf;
% d:      postIve value is delay； negtive is advance；
% 正值表示原始数据有时延，即前面有部分数据不需要，采样点应向后延d再开始取；
% 负值表示原始数据有提前，即前面有部分数据没有取到，采样点应向前提前d再开始取；


function [foe_crs d] = TimingFreqOffset(H_CRS)

global sys_para;

% CFO

H_cor = [H_CRS(:,1) H_CRS(:,3);H_CRS(:,2) H_CRS(:,4)];
xcor_Hls = sum(H_cor(:,2).*conj(H_cor(:,1)));

foe_crs = angle(xcor_Hls)/pi*sys_para.sub_carrier_spacing/15;          %[-1000Hz 1000Hz]


%CTO
m = 48;
H1 = H_CRS(1:200-m/6,:);
H2 = H_CRS(1+m/6:200,:);
xcor_to = sum(sum(H1.*conj(H2),1),2);
d = angle(xcor_Hls)/2/pi/m*sys_para.N_FFT;



