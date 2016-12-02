
function [sbase] = sbaseGen

% nID2 = NID2;
x = zeros(1,31);
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
sbase(i) = 1-2*x(i);



