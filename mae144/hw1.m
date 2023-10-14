%% Homework 1
%% Problem 1

% https://github.com/arrudafar

%% Problem 2a
 
% tf=G*D/(1+G*D)
% G = b/a and D = y/x and so tf= b*y/(a*x + b*y)
% so we get a(s)*x(s)+b(s)*y(s)=f(s)
% Gs = [(s+2)(s-2)(s+5)(s-5)]/[(s+1)(s-1)(s+3)(s-3)(s+6)(s-6)] =b/a
b=RR_poly([2 -2 5 -5],1)
a=RR_poly([1 -1 3 -3 6 -6],1)
% our desired target for system using controller that puts poles at 
% s={-1,-1,-3,-3,-6,-6}
f0=RR_poly([-1 -1 -3 -3 -6 -6],1)
[x0,y0] = RR_diophantine(a,b,f0)
test0=trim(a*x0+b*y0)
residual0=norm(f0-test0)
% by this solution we see the order of the x and y polys are: 
% x = 3 < y = 5 making it so this solution is not proper and its 
% corresponding residual is at 1.13e-13.

%% problem 2b

b=RR_poly([2 -2 5 -5],1)
a=RR_poly([1 -1 3 -3 6 -6],1)
% lets add some s=-20 poles on my f target
f1=RR_poly([-1 -1 -3 -3 -6 -6 -20],1)
[x1,y1] = RR_diophantine(a,b,f1)
test1=trim(a*x1+b*y1)
residual1=norm(f1-test1)
% when we add 1 -20 pole
f2=RR_poly([-1 -1 -3 -3 -6 -6 -20 -20],1)
[x2,y2] = RR_diophantine(a,b,f2)
test2=trim(a*x2+b*y2)
residual2=norm(f2-test2)
% when we add 2 -20 pole
f3=RR_poly([-1 -1 -3 -3 -6 -6 -20 -20 -20],1)
[x3,y3] = RR_diophantine(a,b,f3)
test3=trim(a*x3+b*y3)
residual3=norm(f3-test3)
% when we add 3 -20 pole
f4=RR_poly([-1 -1 -3 -3 -6 -6 -20 -20 -20 -20],1)
[x4,y4] = RR_diophantine(a,b,f4)
test4=trim(a*x4+b*y4)
residual4=norm(f4-test4)
% when we add 6 -20 pole
f5=RR_poly([-1 -1 -3 -3 -6 -6 -20 -20 -20 -20 -20 -20],1)
[x5,y5] = RR_diophantine(a,b,f5)
test5=trim(a*x5+b*y5)
residual5=norm(f5-test5)
% we notice that only when we add about 6 -20 poles to our target that 
% the solution becomes proper with x = 6 > y = 5 with a residual at 
% 2.5e-04 which is still close to zero.

%% problem 3

% to perform a matched z transform we want to replace each first order term
% in the form of [s+a] by its digital equivalent of [1-(e^(-aT))*z^(-1))]
% so the idea is we take a transfer function of the from of:
% ts = [(s+b1)(s+b2)...(s+bm)]/[(s+a1)(s+a2)...(s+an)] =b(s)/a(s)
% and turn it into its matched z transformation.

% here are the tests given for the C2D_zoh function
bs=[1]; as=[1 2 1]; h=0.01; 
Gs=RR_tf(bs,as)
[Gz]=RR_C2D_zoh(Gs,h)
disp('Corresponding Matlab solution:'), c2d(tf(bs,as),h,'zoh')

% here are the tests given for the C2D_tustin function
ys=20*[1 1]; xs=[1 10]; h=0.01; Ds=RR_tf(ys,xs); omegac=sqrt(10); 
[Dz]=RR_C2D_tustin(Ds,h,omegac)
disp('Corresponding Matlab solution:')
opt = c2dOptions('Method','tustin','PrewarpFrequency',omegac); 
c2d(tf(ys,xs),h,opt)
