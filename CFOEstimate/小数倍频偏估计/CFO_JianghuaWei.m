% close all; 
% clear all; 
% clc; 
%参数定义 
function [FFO_df] = CFO_JianghuaWei(N,Ng,recv)
d = Ng+1;
FFO_df1 =  get_FFO_roughly(recv,d,N);
recv_revised = FFO_repair(recv,FFO_df1,N);
FFO_df2 = get_FFO_nicely(recv_revised,d,N);
FFO_df = FFO_df1 + FFO_df2;
%FFO_df1
%FFO_df2
end
%============小数倍频偏粗估计=================
%   FFO_df:小数倍频偏
%   recv:接收到的数据（模拟加噪，频偏后的数据）
%   d:准确的符号定时点
%   N:OFDM符号长度（FFT点数）
function [FFO_df] = get_FFO_roughly(recv,d,N)
P = 0;
for m=0:N/4-1
     P = P + recv(d+m)*conj(recv(d+N/4+m));
end
FFO_df= -2*angle(P)/pi;%小数频偏估计范围为[-1,1]
end
%============小数倍频偏粗估计=================
%   FFO_df:小数倍频偏
%   recv:接收到的数据（模拟加噪，频偏后的数据）
%   d:准确的符号定时点
%   N:OFDM符号长度（FFT点数）
function [FFO_df] = get_FFO_nicely(recv,d,N)
P = 0;
for m=0:N/4-1
     P = P + recv(d+m)*conj(recv(d+3*N/4+m));
end
FFO_df= -2*angle(P)/(3*pi);%小数频偏估计范围为[-1,1]
end
%============小数数倍频偏校正=================
%   FFO_df:小数倍频偏
%   recv:接收到的数据（模拟加噪，频偏后的数据）
%   N:OFDM符号长度（FFT点数）
%   recv_revised:小数倍频偏校正后的数据
function [recv_revised] = FFO_repair(recv,FFO_df,N)
    for k=1:length(recv)
        recv(k) = recv(k)*exp(-sqrt(-1)*2*pi*FFO_df*k/N);
    end
    recv_revised = recv;
end

