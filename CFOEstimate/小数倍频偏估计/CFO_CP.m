% close all; 
% clear all; 
% clc; 
%基于循环前缀进行小数倍频偏估计，该算法只能进行小数倍频偏的估计，不能进行整数倍频偏的估计。
function [FFO_df] = CFO_CP(N,Ng,recv)
d=1;
FFO_df = get_FFO(recv,d,Ng,N);
end
%============频偏估计=====================
function [df] = get_FFO(recv,d,Ng,N)
    P = 0;
    for m=0:Ng-1
        P = P + conj(recv(d+m))*recv(d+m+N);
    end
    df = angle(P)/(2*pi);
end
