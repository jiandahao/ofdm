% close all; 
clear all; 
clc; 
%============获取OFDM系统配置==========
config = OFDMSystemConfig;
N = config.N;
ffo_df = config.ffo_df;
ifo_df = config.ifo_df;
SNR = config.SNR;
simu_times = config.simu_times;

Ng = [N/2 N/4 N/8 N/16];
%=============生成数据符号=============
QAMTable =[7+7i,-7+7i,-7-7i,7-7i]; 
src = QAMTable(randi([0,3],N,1)+1); 
sym = ifft(src,N); 
%cp_sym=sym;%[sym(1,N-Ng+1:N) sym];
recv_sig = zeros(length(Ng),N/2+N);
for n = 1:length(Ng)
    recv_sig(n,1:N+Ng(n)) = [sym(1,N-Ng(n)+1:N) sym];
end
%==============计算Ng = N/8，无噪声情况下频偏估计的情况================
deltaf =  -1.5:0.05:1.5;
CFO_CP_df = zeros(1,length(deltaf));
m = 3;
CP_data = zeros(1,N+Ng(m));
for n=1:length(deltaf)  
    for k=1:N+Ng(m)
        %添加频偏
        CP_data(k) = recv_sig(m,k)*exp(sqrt(-1)*2*pi*deltaf(n)*k/N);
    end
    CFO_CP_df(n) = CFO_CP(N,Ng(m),CP_data);
end
%==============绘制不同频偏值下算法频偏估计值曲线（Ng=N/8,无噪声）====
figure
plot(deltaf,CFO_CP_df,'LineWidth',1);
grid on
set(gca,'XTick',-1.5:0.5:1.5);
xlabel('归一化频偏'); 
ylabel('Van de Beek J J算法频偏估计均值'); 
%================计算不同CP长度，不同信噪比环境下，频偏估计均方误差（deltaf=0.3）=====================
CFO_CP_df_mse = zeros(length(Ng),length(SNR));
CFO_CP_df_mfo = zeros(length(Ng),length(SNR));
%simu_times = 10000;
df = ffo_df+ifo_df;
CP_max_offset = zeros(length(Ng),length(SNR));%记录最大偏差率。
for m=1:length(Ng)
    CP_data = zeros(1,N+Ng(m));
    for k=1:N+Ng(m)
        %添加频偏
        CP_data(k) = recv_sig(m,k)*exp(sqrt(-1)*2*pi*df*k/N);
    end    
    for n=1:length(SNR)  
        for s = 1:simu_times
            CP_data_add_noise = awgn(CP_data,SNR(n));
            df_temp = CFO_CP(N,Ng(m),CP_data_add_noise);
            CFO_CP_df_mse(m,n) = CFO_CP_df_mse(m,n)+(df_temp - df)^2;
            CFO_CP_df_mfo(m,n) = CFO_CP_df_mfo(m,n)+ df_temp;
            CP_max_offset(m,n) = max(CP_max_offset(m,n),abs(df_temp-df));
        end
        CFO_CP_df_mse(m,n) = CFO_CP_df_mse(m,n) / simu_times;
        CFO_CP_df_mfo(m,n) = CFO_CP_df_mfo(m,n) / simu_times;
    end
end
%==========================绘制MSE曲线=========================
figure
d = 1:length(SNR);
label = ['-oc';'-db';'-^k';'-vr'];
for m = 1:length(Ng)
   % plot(SNR,CFO_CP_df_mse(m,d),label(m,:),'LineWidth',1);
   semilogy(SNR,CFO_CP_df_mse(m,d),label(m,:));
    hold on
end
grid on
set(gca,'XTick',0:5:25);

legend('CP=128','CP=64','CP=32','CP=16');
xlabel('SNR(dB)'); 
ylabel('Van de Beek J J算法频偏估计均方误差（MSE）'); 
%==========================绘制MFO曲线========================
% figure
% d = 1:length(SNR);
% label = ['-oc';'-db';'-^k';'-vr'];
% for m = 1:length(Ng)
%     plot(SNR,CFO_CP_df_mfo(m,d),label(m,:));
%     hold on
% end
% grid on
% set(gca,'XTick',0:5:25);
% legend('CP=128','CP=64','CP=32','CP=16');
% xlabel('SNR(dB)'); 
% ylabel('Van de Beek J J算法频偏估计均值（MFO）'); 
%====================将MSE结果记录在文件中================
save('CFO_CP_df_mse.mat','CFO_CP_df_mse','CFO_CP_df_mfo','CP_max_offset');