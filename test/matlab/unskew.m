function [vec]= unskew(a)
vec(1) = a(3,2);
vec(2) = a(1,3);
vec(3) = a(2,1);
end
