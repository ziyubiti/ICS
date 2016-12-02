
% single tone test vector

f = 1e6;
fs = 30.72e6;
ni = 0:fs/1000-1;
t = ni/fs;
cosTV = cos(2*pi*f*t);
sinTV = sin(2*pi*f*t);

tv = cosTV + 1j*sinTV;
plot(abs(tv));grid on;


Fcos = fft(tv)/30720;


figure;plot(abs(Fcos),'-');grid on;


 % TV;Q(1 0 15)
 TV = round(tv*2.^15);
 


Ifile = 'Single_tone_TV.txt';

fid1 = fopen(Ifile,'w+');
for i = 1:1:size(TV,2)
    fprintf(fid1,'%d ,    %d ,\n' ,real(TV(1,i)),imag(TV(1,i)));    
end;






