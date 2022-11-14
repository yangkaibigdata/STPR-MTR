function [f,g]=ob(t)
 x=t(1);
% y=t(2);

%z =t(3);
%h =t(4); 
f = 1/x;
if nargout>1
    g = -(1/(x^2));
end
%w = log(1/(2-x));
%f=(x-3)*(x-3)*source+(y+1)*(y+1)*target+(z-4)*(z-4)+(h-5)*(h-5);
%f = w+y*y*source+z*z+h*h*target;
%f=(x-2)*x*source+x*y*target-(y-1)*y;
end
