function [N,Ng,ffo_df,ifo_df,SNR,QAMTable,simu_times] = get_OFDMSystem_config()
ffo_df = 0.5;%小数倍频偏
ifo_df = 3;%整数倍频偏
simu_times = 1000;%仿真次数
N = 256;%OFDM符号周期长度
Ng = N/8;%循环前缀长度
SNR = [0 5 10 15 20 25];%信噪比
QAMTable =[7+7i,-7+7i,-7-7i,7-7i]; %QAM调制
end