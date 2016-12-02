function correlation_circle = circle_conv(y,x)

x = squeeze(x).';
y = squeeze(y);

if length(x)<length(y)
    z = x;
    x = y;
    y = z;
end
length_x = length(x);
length_y = length(y);
x = [x x(1:length(y)-1)]; 
for i = 1:length_x
    correlation_circle(i) = sum(x(i:length_y+i-1).*conj(y));
end
% output = [output(2:end) output];


