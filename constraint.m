function [c,ceq]=constraint(x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
c1=-x(1);
c2=-x(2);
c3=1-x(3);
c4=-x(4);
c5=-x(5);
c6=-x(6);
c7=-x(7);
c=[c1;c2;c3;c4;c5;c6;c7];%c1];
ceq=[];
end

