
function [c0, c1] = c0c1Gen(NID2)

nID2 = NID2;
x = zeros(1,31);
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


