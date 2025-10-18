
N=200;
x=1:N;
a=0.8;
b=0.3;
c=50;
y=a./(1+exp(-b*(x-c)))+0.1*wgn(1,N,1);

x=x';
y=y';
figure; hold on;
scatter(x,y)

[coeffs,gof]=gng_sigmoid_fit(x,y,[0.75 1 c+5 0.05])

plot(coeffs.a./(1+exp(-coeffs.b*(x-coeffs.c)))+coeffs.d,'r')

%%

y=D14_lc_c_al4(11,:);
%y=D14_lc_nc_al2(26,:);
N=length(y);
inds=~isnan(y);
ind1=find(inds==1,1);
ind2=find(inds==1, 1, 'last' );

a=0.7;
b=0.3;
c=N/2-ind1;
d=0.05;

x=1:(ind2-ind1+1);
y=y(ind1:ind2);

figure; hold on;
scatter(x,y)

[coeffs,gof]=gng_sigmoid_fit(x,y,[a b c+5 d]);
plot(coeffs.a./(1+exp(-coeffs.b*(x-coeffs.c)))+coeffs.d,'r')