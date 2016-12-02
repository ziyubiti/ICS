function output = circle_corr(x,y)
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
% y = conj(y);
length_y = length(y);
for i = 1:length(x)/2
    output(i) = sum(x(i:length_y+i-1).*y);
end
% output = [output(2:end) output];





