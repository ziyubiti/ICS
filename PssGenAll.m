% PSS signal

function [du, fre_sequence, time_sequence] = PssGenAll(n_FFT)


N_FFT = n_FFT;     % 1.4M带宽，下采样16倍,128点
%u = rootNo;  % root index,u=25/29 or 34;

[du1 fre_seq1 time_seq1] = PssGen(N_FFT,25);
[du2 fre_seq2 time_seq2] = PssGen(N_FFT,29);
[du3 fre_seq3 time_seq3] = PssGen(N_FFT,34);

du = [du1.'; du2.'; du3.'];
fre_sequence = [fre_seq1.'; fre_seq2.'; fre_seq3.'];
time_sequence = [time_seq1.'; time_seq2.'; time_seq3.'];

fre2 = sqrt(1.0/N_FFT)*fft(time_sequence,64,2);

% 
% % % A为频域相关矩阵
% A{1,2} = abs(circle_corr(du(1,:),du(2,:)));
% figure;stem(A{1,2});
% A{1,3} = abs(circle_corr(du(1,:),du(3,:)));
% figure;stem(A{1,3});
% A{2,3} = abs(circle_corr(du(2,:),du(3,:)));
% figure;stem(A{2,3});
% A{1,1} = abs(circle_corr(du(1,:),du(1,:)));
% figure;stem(A{1,1});
% A{2,2} = abs(circle_corr(du(2,:),du(2,:)));
% figure;stem(A{2,2});
% A{3,3} = abs(circle_corr(du(3,:),du(3,:)));
% figure;stem(A{3,3});
% % % 
% A{1,2} = abs(circle_corr(fre_sequence(1,:),fre_sequence(2,:)));
% figure;stem(A{1,2});
% A{1,3} = abs(circle_corr(fre_sequence(1,:),fre_sequence(3,:)));
% figure;stem(A{1,3});
% A{2,3} = abs(circle_corr(fre_sequence(2,:),fre_sequence(3,:)));
% figure;stem(A{2,3});
% A{1,1} = abs(circle_corr(fre_sequence(1,:),fre_sequence(1,:)));
% figure;stem(A{1,1});
% A{2,2} = abs(circle_corr(fre_sequence(2,:),fre_sequence(2,:)));
% figure;stem(A{2,2});
% A{3,3} = abs(circle_corr(fre_sequence(3,:),fre_sequence(3,:)));
% figure;stem(A{3,3});
% % % 
% % % 
% % % 
% % % B为时域相关矩阵
% B{1,2} = abs(circle_corr(time_sequence(1,:),time_sequence(2,:)));
% figure;stem(B{1,2});
% B{1,3} = abs(circle_corr(time_sequence(1,:),time_sequence(3,:)));
% figure;stem(B{1,3});
% B{2,3} = abs(circle_corr(time_sequence(2,:),time_sequence(3,:)));
% figure;stem(B{2,3});
% B{1,1} = abs(circle_corr(time_sequence(1,:),time_sequence(1,:)));
% figure;stem(B{1,1});
% B{2,2} = abs(circle_corr(time_sequence(2,:),time_sequence(2,:)));
% figure;stem(B{2,2});
% B{3,3} = abs(circle_corr(time_sequence(3,:),time_sequence(3,:)));
% figure;stem(B{3,3});
% 
% % 


