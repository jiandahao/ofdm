%--carrier frequency recovery-- 
%--MM algorithm-- 
%--2010/09/06-- 
clc 
clear all 
close all 
size_pilot = 36; 
num_pilot = 50; 
freq_offset =0.05;       %归一化频偏 
phase_offset =0;  %rand(1)*2*pi; 
SNR = 0; 
temp_cap = 0; 
auto_correlation = zeros(1,18); 
freq_cap = zeros(1,num_pilot);%%用来计算每个导频块所捕获的频偏(归一化频偏) 
ca0 = freq_cap; 
signal = ones(1,num_pilot*size_pilot)*exp(j*pi/4); 
signal_chan = channel(SNR, signal, freq_offset, phase_offset); 
%  plot(signal_chan,'*');   
ii = 1; 
flag = 0; 
while ii <= num_pilot*size_pilot 
    flag = flag+1; 
    %%--------MM算法---------- 
    data = signal_chan(ii:ii+35).*exp(j*pi*7/4); 
    %%计算自相关函数 
    for m = 1:size_pilot/2 
        for k= (m+1):length(data) 
            auto_correlation(m) = auto_correlation(m) + exp(j*(angle(data(k))-angle(data(k-m)))); 
        end 
        auto_correlation(m) = auto_correlation(m)/(length(data)-m); 
    end 
    %%计算平滑校正系数 
    for m = 1:size_pilot/2 
        coef(m) = (3*((36-m)*(36-m+1)-18*(36-18)))/(18*(4*18^2-6*18*36+3*36^2-1)); 
    end 
    %%计算频偏 
    for m = 2:size_pilot/2 
        temp = angle(auto_correlation(m))-angle(auto_correlation(m-1)); 
        freq_cap(flag) = freq_cap(flag) + coef(m)*temp; 
    end 
    freq_cap(flag)=freq_cap(flag)/(2*pi); %归一化频偏   
    ii = ii+ 36; 
    for jj =ii:num_pilot*size_pilot 
        signal_chan(jj) = signal_chan(jj)*exp(-j*2*pi*(freq_cap(flag)).*(jj-1)); 
    end 
end 
signal_out = signal_chan; 
% signal_out = angle(signal_chan)/(2*pi); 
freq_cap; 
figure(1) 
plot(signal_out,'*'); 
title ('受到频偏干扰的曲线'); 
figure(2) 
plot(freq_cap(1:num_pilot)); 
grid on; 
title ('捕获的频偏'); 
% axis([1 num_pilot -0.2 0.2]); 
% axis([1 36*num_pilot -2 2]); 