%该类用于设置OFDM系统参数，用于仿真使用
classdef OFDMSystemConfig
    %Don't modify the parameter's name of the following properties.
   properties
    ffo_df =0.4;%小数倍频偏
    ifo_df =16;%整数倍频偏
    simu_times = 1000;%仿真次数
    N = 256;%OFDM符号周期长度
    Ng = 32;%循环前缀长度,N/8
    SNR = [-4 -2 0 2 4 6 8 10 12 14 16 20];%信噪比
    %SNR = 0;
    %QAMTable =[7+7i,-7+7i,-7-7i,7-7i]; %QAM调制
   end 
%    methods
%        function obj = OFDMSystemConfig()
%            obj.ffo_df =rand(1,obj.simu_times)-0.5;%小数倍频偏
%        end
%    end
end
