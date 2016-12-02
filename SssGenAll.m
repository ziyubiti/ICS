
function [freq_sss_all, time_sss_all] = SssGenAll(n_FFT,NID2)

N_FFT = n_FFT;
% NID2 = 0;
% NID1 = 0;
freq_sss_all = zeros(336,62);  %168*2 =336, first 168 line are sf0 data; second 168 line are sf 5 data;
time_sss_all = zeros(336,64);

for NID1 = 0:1:167
    [freq_sss, freq_sss_even, freq_sss_odd,time_sss] = SssGen(N_FFT,NID2,NID1);
    freq_sss_all(NID1+1,:) = freq_sss(1,:);
    freq_sss_all(NID1+1+168,:) = freq_sss(2,:);
    time_sss_all(NID1+1,:) = time_sss(1,:);
    time_sss_all(NID1+1+168,:) = time_sss(2,:);
end

%% SSS para calc
% n = 0:167;
% qh = floor(n/30);
% q = floor((n + qh.*(qh+1)/2)/30);
% mh = n + q.*(1+q)/2;
% m0 = mod(mh,31);
% m1 = mod((m0+floor(mh./31)+1),31);
% SSSpara = [n.' qh.' q.' mh.' m0.' m1.'];
% 
% chat = [ 0,0,0,0,1,0,1,0,1,1,1,0,1,1,0,0,0,1,1,1,1,1,0,0,1,1,0,1,0,0,1];
% ch = 1 - 2 * chat;
% chloop = [ch ch];
% shat = [0,0,0,0,1,0,0,1,0,1,1,0,0,1,1,1,1,1,0,0,0,1,1,0,1,1,1,0,1,0,1];
% sh = 1 - 2 * shat;
% shloop = [sh sh];
% zhat = [0,0,0,0,1,1,1,0,0,1,1,0,1,1,1,1,1,0,1,0,0,0,1,0,0,1,0,1,0,1,1];
% zh = 1 - 2 * zhat;
% zhloop = [zh zh];




