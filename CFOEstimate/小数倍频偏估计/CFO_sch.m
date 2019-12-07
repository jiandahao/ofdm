% close all; 
% clear all; 
% clc; 
%参数定义 
function [FFO_df,IFO_df,CFO_df] = CFO_sch(N,Ng,recv,v_sig)
d = Ng+1;
FFO_df =  get_FFO(recv,d,N);
recv_revised = FFO_repair(recv,FFO_df,N);
IFO_df = get_IFO(recv_revised,d,N,Ng,v_sig);
CFO_df = 2*IFO_df + FFO_df;
end
%============小数倍频偏估计=================
%   FFO_df:小数倍频偏
%   recv:接收到的数据（模拟加噪，频偏后的数据）
%   d:准确的符号定时点
%   N:OFDM符号长度（FFT点数）
function [FFO_df] = get_FFO(recv,d,N)
P = 0;
for m=0:N/2-1
     P = P + conj(recv(d+m))*recv(d+N/2+m);
end
FFO_df= angle(P)/pi;%小数频偏估计范围为[-1,1]
end

%============小数数倍频偏校正=================
%   FFO_df:小数倍频偏
%   recv:接收到的数据（模拟加噪，频偏后的数据）
%   N:OFDM符号长度（FFT点数）
%   recv_revised:小数倍频偏校正后的数据
function [recv_revised] = FFO_repair(recv,FFO_df,N)
    for k=1:length(recv)
        recv(k) = recv(k)*exp(-i*2*pi*FFO_df*k/N);
    end
    recv_revised = recv;
end
%============整数数倍频偏估计=================
%   N:OFDM符号长度（FFT点数）
%   Ng:循环前缀长度
%   recv_revised:小数倍频偏校正后的数据
%   d:准确的符号定时点
%   v:差分序列
%   FFO_df:整数倍频偏估计值
function [IFO_df] =  get_IFO(recv_revised,d,N,Ng,v)
recv = recv_revised;
Ns=N+Ng;
r1 = recv(1,Ng+1:Ns);%获取第一个OFDM符号
r2 = recv(1,N+2*Ng+1:2*Ns);%获取第二个OFDM符号
R1 = fft(r1);
R2 = fft(r2);

for g = -10:10
    for k = 1:2:N
        if k+2*g <=0
             B1(k) = conj(R1(k+2*g+N))*conj(v(k))*R2(k+2*g+N);
        end
        if (k+2*g<=N) &&  (k+2*g>0)
            B1(k) = conj(R1(k+2*g))*conj(v(k))*R2(k+2*g);
        end
        if k+2*g>N
            B1(k) = 0;
        end
        B2(k) = abs(R2(k))^2;
    end
    B(g+11) = (abs(sum(B1)))^2/(2*sum(B2)^2);
end
[~,IFO_df] = max(B);
IFO_df = IFO_df - 11;
%CFO_df = 2*(IFO_df-11) + FFO_df
end

% R1 = fft(r1);
% R2 = fft(r2);
% P = zeros(1,N);
% R = zeros(1,N);
% B =zeros(1,10);
% for g=0:9
%     for k=1:2:N
%         if k+2*g<=N
%              P(g+1) =  P(g+1)+ conj(R1(k+2*g))*conj(v(k))*R2(k+2*g);
%         end
%         if k+2*g>N
%             P(g+1)=P(g+1)+0;
%         end
%         R(g+1) =  R(g+1)+abs(R2(g+1))^2;
%     end
%     B(g+1) = power(abs(P(g+1)),2)/(2*power(R(g+1),2));
% end
% d=1:10;
% plot(d,B)

