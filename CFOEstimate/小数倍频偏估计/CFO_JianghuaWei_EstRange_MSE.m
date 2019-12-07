clear all; 
clc; 
%============获取OFDM系统配置==========
config = OFDMSystemConfig;
N = config.N;
Ng = config.Ng;
ffo_df = config.ffo_df;
ifo_df = config.ifo_df;
SNR = config.SNR;
simu_times = config.simu_times;
%=============生成数据符号=============
QAMTable =[7+7i,-7+7i,-7-7i,7-7i]; 
src = QAMTable(randi([0,3],N,1)+1); 
sym = ifft(src,N); 
cp_sym = [sym(1,N-Ng+1:N) sym];
%=============生成训练符号==============
buf=QAMTable(randi([0,3],N/4,1)+1); %1x(N/2)的矩阵
train1=zeros(1,N/2); 
index = 1; 
for n=1:2:N/2 
     train1(n)=buf(index); 
     index=index+1; 
end; 
sch = ifft(train1,N/2);   %[A A]的形式 
cp_train = [sch(1,N/4-Ng+1:N/4) sch sch];%[cp A A A A]格式
transmit_data = [cp_train cp_sym];
%==============计算无噪声情况下频偏估计的情况================
deltaf = -3:0.05:3;
CFO_Jiang_df = zeros(1,length(deltaf));
recv_sig = zeros(1,length(transmit_data));
for n=1:length(deltaf)  
    for k=1:length(transmit_data)
        %添加频偏
        recv_sig(k) = transmit_data(k)*exp(sqrt(-1)*2*pi*deltaf(n)*k/N);
    end
    CFO_Jiang_df(n) =  CFO_JianghuaWei(N,Ng,recv_sig);
end
%==============绘制不同频偏值下算法频偏估计值曲线（无噪声）====
figure
plot(deltaf,CFO_Jiang_df,'LineWidth',1);
grid on
set(gca,'XTick',-3:0.5:3);
xlabel('归一化频偏'); 
ylabel('Moose算法频偏估计均值'); 
%================计算不同信噪比环境下，频偏估计均方误差（deltaf=0.3）=====================
CFO_Jiang_df_mse = zeros(1,length(SNR));
CFO_Jiang_df_mfo = zeros(1,length(SNR));
%simu_times = 10000;
df = ffo_df+ifo_df;
recv_sig = zeros(1,length(transmit_data));
Jiang_max_offset = zeros(1,length(SNR));%记录最大偏差率。
for k=1:length(transmit_data)
    %添加频偏
    recv_sig(k) = transmit_data(k)*exp(sqrt(-1)*2*pi*df*k/N);
end    
for n=1:length(SNR)  
    for s = 1:simu_times
        recv_sig_add_noise = awgn(recv_sig,SNR(n));
        df_temp = CFO_JianghuaWei(N,Ng,recv_sig_add_noise);
        CFO_Jiang_df_mse(n) = CFO_Jiang_df_mse(n)+(df_temp - df)^2;
        CFO_Jiang_df_mfo(n) = CFO_Jiang_df_mfo(n)+ df_temp;
        Jiang_max_offset(n) = max(Jiang_max_offset(n),abs(df_temp-df));
    end
    CFO_Jiang_df_mse(n) = CFO_Jiang_df_mse(n) / simu_times;
    CFO_Jiang_df_mfo(n) = CFO_Jiang_df_mfo(n) / simu_times;
end
%==========================绘制MSE曲线=========================
figure
d = 1:length(SNR);
semilogy(SNR,CFO_Jiang_df_mse(d));
grid on
set(gca,'XTick',0:5:25);
xlabel('SNR(dB)'); 
ylabel('JianghuaWei算法频偏估计均方误差（MSE）'); 
%====================将MSE结果记录在文件中================
save('CFO_Jiang_df_mse.mat','CFO_Jiang_df_mse','CFO_Jiang_df_mfo','Jiang_max_offset');