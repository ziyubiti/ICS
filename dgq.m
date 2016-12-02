


% dgq tv


dgq = load('C:\receive_file\ltedl20mhziqdata.c');
% 
% for j = 1:2
%     for i = 1:size(dgq,1)
%     dgg_scale2 = dgq*2;
%     dgg_scale4 = dgq*4;
%     end;
% end;
   

dgg_scale2 = dgq*2;
Ifile = 'dgg_scale2.txt';

fid1 = fopen(Ifile,'w+');
for i = 1:1:size(dgg_scale2,1)
    fprintf(fid1,'%d     %d \n' ,dgg_scale2(i,1),dgg_scale2(i,2));    
end;

