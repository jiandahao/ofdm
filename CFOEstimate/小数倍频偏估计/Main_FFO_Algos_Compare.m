close all;
clear all;
clc

% CFO_sch_EstRange_MSE;
% CFO_moose_EstRange_MSE;
% CFO_CP_EstRange_MSE;
% CFO_JianghuaWei_EstRange_MSE;
% CFO_proposed_EstRange_MSE;
clear all;
clc
config = OFDMSystemConfig;
SNR = config.SNR;
d = 1:length(SNR);
load CFO_Jiang_df_mse.mat;
load CFO_sch_df_mse.mat;
load CFO_moose_df_mse.mat;
load CFO_CP_df_mse.mat
load CFO_proposed_df_mse.mat;
%=============绘制最大偏差率曲线（Maximum offset rate）========
figure
plot(SNR,CP_max_offset(3,d),'-oc','LineWidth',1);
hold on
plot(SNR,sch_max_offset(d),'-db','LineWidth',1);
hold on
plot(SNR,moose_max_offset(d),'-^k','LineWidth',1);
hold on
plot(SNR,Jiang_max_offset(d),'-*','LineWidth',1);
hold on
plot(SNR,proposed_max_offset(d),'-vr','LineWidth',1);
hold on
grid on
set(gca,'XTick',SNR);
set(gca,'YTick',0:0.1:1.6);
xlabel('SNR(dB)'); 
ylabel('最大估计偏差值'); 
legend('Van de Beek算法(Ng=32,N=256)','Schmidl&Cox算法','Moose算法','Jianghua Wei算法','本文算法');
grid on;
%=========绘制CP算法、moose算法与Schmidl算法MSE曲线========
figure 
semilogy(SNR,CFO_CP_df_mse(3,d),'-vr','LineWidth',1);
hold on
semilogy(SNR,CFO_moose_df_mse(d),'-db','LineWidth',1);
hold on
semilogy(SNR,CFO_sch_df_mse(d),'-^k','LineWidth',1);
set(gca,'XTick',SNR);
xlabel('SNR(dB)'); 
ylabel('频偏估计均方误差（MSE）'); 
legend('Van de Beek算法(Ng=32,N=256)','Moose算法','Schmidl&Cox算法');
grid on;

%=========绘制Schmidl算法与Jianghua Wei算法MSE曲线对比图======
figure 
semilogy(SNR,CFO_sch_df_mse(d),'-^k','LineWidth',1);
hold on
semilogy(SNR,CFO_Jiang_df_mse(d),'-vr','LineWidth',1);
hold on
set(gca,'XTick',SNR);
xlabel('SNR(dB)'); 
ylabel('频偏估计均方误差（MSE）'); 
legend('Schmidl&Cox算法','Jianghua Wei算法');
grid on;
%===========绘制频偏估计均方误差曲线==============
figure
semilogy(SNR,CFO_CP_df_mse(3,d),'-oc','LineWidth',1);
hold on
semilogy(SNR,CFO_sch_df_mse(d),'-db','LineWidth',1);
hold on
semilogy(SNR,CFO_moose_df_mse(d),'-^k','LineWidth',1);
hold on
semilogy(SNR,CFO_Jiang_df_mse(d),'-*','LineWidth',1);
hold on
semilogy(SNR,CFO_proposed_df_mse(d),'-vr','LineWidth',1);
hold on
grid on
set(gca,'XTick',SNR);
xlabel('SNR(dB)'); 
ylabel('频偏估计均方误差（MSE）'); 
legend('Van de Beek算法(Ng=32,N=256)','Schmidl&Cox算法','Moose算法','Jianghua Wei算法','本文算法');
%===========绘制频偏估计均值曲线==============
figure
plot(SNR,CFO_CP_df_mfo(3,d),'-oc','LineWidth',1);
hold on
plot(SNR,CFO_sch_df_mfo(d),'-db','LineWidth',1);
hold on
plot(SNR,CFO_moose_df_mfo(d),'-^k','LineWidth',1);
hold on
plot(SNR,CFO_Jiang_df_mfo(d),'-*','LineWidth',1);
hold on
plot(SNR,CFO_proposed_df_mfo(d),'-vr','LineWidth',1);
hold on
grid on
set(gca,'XTick',SNR);
xlabel('SNR(dB)'); 
ylabel('小数倍频偏估计均值'); 
legend('Van de Beek算法(Ng=32,N=256)','Schmidl&Cox算法','Moose算法','Jianghua Wei算法','本文算法');
% SNR = 0:2:24;
% d = 1:length(SNR);
% figure
% semilogy(SNR,CFO_moose_df_mse(d),'-vr','LineWidth',1);
% hold on
% semilogy(SNR,CFO_CP_df_mse(3,d),'-^k','LineWidth',1);
% hold on;
% semilogy(SNR,CFO_sch_df_mse(d),'-x','LineWidth',1);
% figure
% plot(SNR,moose_df(d),'-vr','LineWidth',1);
% hold on
% plot(SNR,CP_FFO_df(d),'-oc','LineWidth',1);
% hold on
% plot(SNR,sch_FFO_df(d),'-x','LineWidth',1);
% hold on
% plot(SNR,CFO_Jiang_df_mfo(d),'-*','LineWidth',1);
% hold on
% plot(SNR,p_FFO_df(d),'-^k','LineWidth',1);
% 
% hold on
% xlabel('SNR(dB)'); 
% ylabel('频偏估计均方误差（MSE）'); 
% legend('Moose算法','Van de Beek JJ算法','Schmidl&Cox算法','Jianghua Wei算法','本文算法');
% grid on;