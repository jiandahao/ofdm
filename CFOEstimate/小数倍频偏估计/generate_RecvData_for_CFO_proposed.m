%函数功能：生成改进算法对应的同步帧序列（TS_sym），以及仿真使用的加入频偏的传输数据（recv_sig）
%输入参数 N:OFDM系统有效符号长度，即FFT电视
%         Ng:循环前缀的符号长度
function [recv_sig,TS_sym] =  generate_RecvData_for_CFO_proposed(N,Ng,deltaf)
%N     %FFT/IFFT 变换的点数或者子载波个数（Nu=N） 
%Ng      %循环前缀的长度 (保护间隔的长度) 
Ns=Ng+N;     %包括循环前缀的符号长度 
QAMTable =[7+7i,-7+7i,-7-7i,7-7i]; 
%=============生成数据符号=============
src = QAMTable(randi([0,3],N,1)+1); 
sym = ifft(src); 
cp_sym=[sym(1,N-Ng+1:N) sym];
predata = cp_sym;
suffixdata = cp_sym;
%=============生成训练序列======================
a = zeros(1,N/4);
mu = N/4-1;
for n=0:N/4-1
   a(n+1) =7*exp(4*sqrt(-1)*mu*pi*n*n/N);
end
A = ifft(a);
%A =a;
B = conj(A(1,N/4:-1:1));
C= zeros(1,N/4);
for n=1:1:N/4 % C(n) = m(n)B(n)
    if mod(n,2)
       % C(n) = (-1)*conj(B(n));
       C(n) = (-1)*B(n);
    else
       % C(n) = conj(B(n));
       C(n) = B(n);
    end
end
D=C;
signal = [A B C D];
TS_sym = signal;
cp_train = [signal(1,N-Ng+1:N) signal];
transmit_data = [predata cp_train suffixdata];
%=============模拟经过信道=========================
recv = transmit_data;
for k=1:length(recv)
    %添加频偏
    recv(k) = recv(k)*exp(i*2*pi*deltaf*k/N);
end
recv_sig = recv;
end