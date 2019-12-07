%函数功能：生成基于循环前缀的频偏估计算法对应仿真使用的加入频偏的传输数据（recv_sig）
%输入参数 N:OFDM系统有效符号长度，即FFT电视
%         Ng:循环前缀的符号长度
function [recv_sig] = generate_RecvData_for_CFO_CP(N,Ng,deltaf)
Ns=Ng+N;     %包括循环前缀的符号长度 
QAMTable =[7+7i,-7+7i,-7-7i,7-7i]; 
%=============生成数据符号=============
src = QAMTable(randi([0,3],N,1)+1); 
sym = ifft(src,N); 
%cp_sym=sym;%[sym(1,N-Ng+1:N) sym];
signal = [sym(1,N-Ng+1:N) sym];
%=============模拟添加频偏和噪声========
recv_sig = signal;
for k=1:length(recv_sig)
    %添加频偏
    recv_sig(k) = recv_sig(k)*exp(sqrt(-1)*2*pi*deltaf*k/N);
end
end