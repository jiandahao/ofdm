% close all; 
% clear all; 
% clc; 
%使用Moose算法进行小数倍频偏估计
function [FFO_df] = CFO_moose(N,Ng,recv)
    d = 0;
    FFO_df  =get_FFO(recv,d,N,Ng);
end
%============频偏估计=====================
function [FFO_df]  = get_FFO(recv,d,N,Ng)
Ns= N+Ng;
r1 = recv(1,Ng+1:Ns);%获取第一个OFDM符号
r2 = recv(1,N+Ng+1:2*N+Ng);%获取第二个OFDM符号
R1 = fft(r1,N);
R2 = fft(r2,N);
P = 0;
for k=1:N
%     Im_sum = Im_sum + imag(R2(k))*conj(imag(R1(k)));
%     Re_sum = Re_sum + real(R2(k))*conj(real(R1(k)));
    P = P + conj(R1(k))*R2(k);
end
FFO_df = angle(P)/(2*pi);%计算频偏，频偏估计范围为[-0.5,0.5]
end


