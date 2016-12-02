
function [z1mi] = z1miGen(mi)

% nID2 = NID2;
x = zeros(1,31);
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
z1mi(n) = z(mod((n-1+ mod(mi,8)),31)+1);


