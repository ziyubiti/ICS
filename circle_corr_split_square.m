
function output = circle_corr_split_square(x,y,M)
% description:
% 计算两个输入数值的循环相关值
% input:
%           x: 1*Nx 
%           y: 1*Ny  (Nx>=Ny)
%              
% output:  
%           output:1*Ny

x = [x x];
x = conj(x);
length_y = length(y);
if (M == 2)
    for i = 1:length(x)/2
        output(i) = abs(sum(x(i:length_y/2+i-1).*y(1,1:length_y/2)))*abs(sum(x(i:length_y/2+i-1).*y(1,1:length_y/2))) ...
            + abs(sum(x(i+length_y/2:length_y+i-1).*y(1,length_y/2+1:length_y)))*abs(sum(x(i+length_y/2:length_y+i-1).*y(1,length_y/2+1:length_y)));
    end;
end;
if (M == 4)
    for i = 1:length(x)/2
        output(i) = abs(sum(x(i:length_y/4+i-1).*y(1,1:length_y/4)))*abs(sum(x(i:length_y/4+i-1).*y(1,1:length_y/4))) ...
            + abs(sum(x(i+length_y/4:length_y/2+i-1).*y(1,length_y/4+1:length_y/2)))*abs(sum(x(i+length_y/4:length_y/2+i-1).*y(1,length_y/4+1:length_y/2)))...
            + abs(sum(x(i+length_y/2:length_y/4*3+i-1).*y(1,length_y/2+1:length_y/4*3)))*abs(sum(x(i+length_y/2:length_y/4*3+i-1).*y(1,length_y/2+1:length_y/4*3)))...
            + abs(sum(x(i+length_y/4*3:length_y+i-1).*y(1,length_y/4*3+1:length_y)))*abs(sum(x(i+length_y/4*3:length_y+i-1).*y(1,length_y/4*3+1:length_y)));
    end;
end;
% output = [output(2:end) output];






