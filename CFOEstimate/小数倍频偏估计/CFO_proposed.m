% close all; 
% clear all; 
% clc; 
%参数定义 
function [FFO_df,IFO_df,CFO_df] =CFO_proposed(N,Ng,recv,Ts_sym)
%添加高斯白噪声
    Ns=N+Ng;
    d = Ns+Ng+1;
    [revised_recv,FFO_df] = FFO_estimate(recv,d,N);
    IFO_df= get_IFO(revised_recv,d,N,Ts_sym);
    
    CFO_df = FFO_df + IFO_df;
end
%============时域频偏估计与修正=====================
function [revised_recv,df] = FFO_estimate(recv,d,N)
%simulation_times = 1;
df1 = 0;
df2 = 0;
%for i=1:simulation_times
    df1 = get_FFO_roughly(recv,d,N);
    recv = FFO_repair(recv,df1,N);
    df2 = get_FFO_finely(recv,d,N);
    recv = FFO_repair(recv,df2,N);
    df = df1+df2;
%end
revised_recv = recv;
end
%===========时域小数倍频偏补偿========================
function [revise_recv] = FFO_repair(recv,df,N)
    for k=1:length(recv)
        recv(k) = recv(k)*exp(-i*2*pi*df*k/N);
    end
    revise_recv = recv;
end
%============时域小数频偏粗估计=====================
%   df:小数倍频偏粗估计值
%   recv:接收到的数据（模拟加噪，频偏后的数据）
%   d:准确的符号定时点
%   N:OFDM符号长度（FFT点数）
function [df] = get_FFO_roughly(recv,d,N)
P = 0;
for m=0:N/4-1
     P = P + conj(recv(d+m+N/2))*recv(d+N/2+N/4+m);
end
df = 2*angle(P)/pi;
end
%============时域小数频偏细估计=====================
%   df:小数倍频偏细估计值
%   recv:小数倍频偏估计值校正后的数据
%   d:准确的符号定时点
%   N:OFDM符号长度（FFT点数）
function [df] = get_FFO_finely(recv,d,N)
P = 0;
for m=0:N/4-1
     P = P + conj(recv(d+m+N/4))*power(-1,m+1)*recv(d+3*N/4+m);
end
df = angle(P)/pi;
end

%============时域整数频偏估计=====================
%参考OFDM同步技术研究及实现——柴立凯
%   df:整数倍频偏估计值
%   FFO_revised_recv:小数倍频偏估计值校正后的数据
%   d:准确的符号定时点
%   N:OFDM符号长度（FFT点数）
%   TS:时域上使用的训练序列
function [df] = get_IFO(FFO_revised_recv,d,N,TS)
P = zeros(1,N/2+1);
for m=-N/4:N/4
    for n=0:N-1
     P(m+N/4+1) = P(m+N/4+1) + FFO_revised_recv(d+n)*conj(TS(n+1))*exp(-sqrt(-1)*2*pi*m*n/N);
    end
end
[~,df] = max(abs(P));
df = df - N/4 - 1;
end