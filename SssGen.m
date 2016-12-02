% SSS signal

% sf 0/5 data,two line;

function [freq_sss, freq_sss_even, freq_sss_odd,time_sss] = SssGen(n_FFT,NID2,NID1)

nID2 = NID2;
nID1 = NID1;
N_FFT = n_FFT;

x = zeros(1,31);
c = zeros(1,31); 
c0 = zeros(1,31); 
c1 = zeros(1,31); 
s = zeros(1,31); 
s0 = zeros(1,31); 
s1 = zeros(1,31); 
z = zeros(1,31); 
z1m0 = zeros(1,31); 
z1m1 = zeros(1,31); 

% step 0: calc m0 and m1; related to NID1;
qt = floor(nID1/30);
q = floor((nID1 + qt*(qt+1)/2)/30);
mt = nID1 + q*(q+1)/2;
m0 = mod(mt,31);
m1 = mod(m0 + floor(mt/31)+1,31);
%     table(nID1+1,1) = nID1;        % for check
%     table(nID1+1,2) = m0;
%     table(nID1+1,3) = m1;



% step 1: c sequence,c0 and c1, related to NID2;
x(1)=0;
x(2)=0;
x(3)=0;
x(4)=0;
x(5)=1;

for i = 1:26
    x(i+5) = mod((x(i+3)+x(i)),2);
end
i = 1:31;
c(i) = 1-2*x(i);
n = 1:31;
c0(n) = c(mod((n-1+ nID2),31)+1);
c1(n) = c(mod((n-1+ nID2 +3),31)+1);

% step 2: s sequence,s0 and s1, related to m0 and m1;
x(1)=0;
x(2)=0;
x(3)=0;
x(4)=0;
x(5)=1;

for i = 1:26
    x(i+5) = mod((x(i+2)+x(i)),2);
end
i = 1:31;
s(i) = 1-2*x(i);
n = 1:31;
s0(n) = s(mod((n-1+ m0),31)+1);
s1(n) = s(mod((n-1+ m1),31)+1);

% step 3: z sequence,z0 and z1, related to m0 and m1;
x(1)=0;
x(2)=0;
x(3)=0;
x(4)=0;
x(5)=1;
for i = 1:26
    x(i+5) = mod((x(i+4)+x(i+2)+x(i+1)+x(i)),2);
end
i = 1:31;
z(i) = 1-2*x(i);
n = 1:31;
z1m0(n) = z(mod((n-1+ mod(m0,8)),31)+1);
z1m1(n) = z(mod((n-1+ mod(m1,8)),31)+1);


%%
n = 1:31;
d(1,2*n-1) = s0.*c0;     % subframe 0;
d(2,2*n-1) = s1.*c0;     % subframe 5;
d(1,2*n) = s1.*c1.*z1m0;     % subframe 0;
d(2,2*n) = s0.*c1.*z1m1;     % subframe 5;

d_even = d(:,2*n-1);
d_odd = d(:,2*n);
freq_sss = d;
freq_sss_even = d_even;
freq_sss_odd = d_odd;

fre_seq_padding = zeros(N_FFT,2);
time_seq = zeros(N_FFT,2);
fre_seq_padding(2:32,:) = d(:,32:62).';    %  µÚ1¸öÎªDC
fre_seq_padding(-30+N_FFT:N_FFT,:) = d(:,1:31).';
time_seq =  sqrt(1.0*N_FFT)*ifft(fre_seq_padding,[],1);   
time_sss = time_seq.';













