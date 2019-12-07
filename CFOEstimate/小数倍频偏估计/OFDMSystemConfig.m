%该类用于设置OFDM系统参数，用于仿真使用
classdef OFDMSystemConfig
    %Don't modify the parameter's name of the following properties.
   properties
    ffo_df = 0.3;%小数倍频偏
    ifo_df = 0;%整数倍频偏
    simu_times = 10000;%仿真次数
    N = 256;%OFDM符号周期长度
    Ng = 32;%循环前缀长度,N/8
    SNR = [0 5 10 15 20 25];%信噪比
     %SNR = 0:2:24;
    %QAMTable =[7+7i,-7+7i,-7-7i,7-7i]; %QAM调制
   end 
end