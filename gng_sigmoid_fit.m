

function [coeffs,gof]=io_sigmoid_fit(x,y,initcond)
%  fittype('a/(1+exp(-b*(x-c)+d)', 'StartPoint', [0.65 1 n 0.05] );

if size(x,1)<size(x,2)
    x=x';
end

if size(y,1)<size(y,2)
    y=y';
end

fo=fitoptions('Method','NonlinearLeastSquares','StartPoint',initcond);
ft=fittype('a/(1+exp(-b*(x-c)))+d','options',fo);
[coeffs,gof]=fit(x,y,ft,'lower',0,'upper',1);
%[coeffs,gof]=fit(x,y,ft);