%close all
clear all
clc
%ggg
%Get OFDM simulation Configs
config = OFDMSystemConfig;
N = config.N;
Ng = config.Ng;
Ns = N+Ng;
ifo_df = config.ifo_df;
simu_times = config.simu_times;
SNR = config.SNR;%
%Initialize Algorithms Class Object
pAlgo = ProposedAlgo;
schAlgo = SchmidlAlgo;
SoumAlgo = SoumitraAlgo;
LiuAlgo = LiuChunAlgo;
RAlgo = RenAlgo;
FangAlgo = FangYiboAlgo;
KimAlgo = YHKimAlgo;
ZhangAlgo = ZhangjieAlgo;
%=====Generate simulation Data for each Algorithms=========
[pData,pTS] = pAlgo.generate_data(N,Ng);
[schData,Sch_v] = schAlgo.generate_data(N,Ng);
[SoumTData,SoumFData] = SoumAlgo.generate_data(N,Ng);
[LiuData,LiuTS] = LiuAlgo.generate_data(N,Ng);
[RenData,RenTS] = RAlgo.generate_data(N,Ng);
[FangData,Fang_v,Fang_mu] = FangAlgo.generate_data(N,Ng);
[KimData,KimDiffSeq] = KimAlgo.generate_data(N,Ng);
[ZhangData,ZhangDiffSeq] = ZhangAlgo.generate_data(N,Ng);
%Arrays for recording Properties of Success (POS) for each Algorithms
pSR = zeros(1,length(SNR));
schSR = zeros(1,length(SNR));
soumSR = zeros(1,length(SNR));
liuSR = zeros(1,length(SNR));
fangSR = zeros(1,length(SNR));
kimSR = zeros(1,length(SNR));
renSR = zeros(1,length(SNR));
zhangSR = zeros(1,length(SNR));
ifo_df = 16;
parfor n = 1:length(SNR)
    for s=1:simu_times
        ffo_df = -0.5 + rand();
        %ffo_df = -0.05 + 0.1*rand();
        %ifo_df = 4;
        df = ifo_df + ffo_df;
        
        %Adding Carrier Frequency Offset into simulation data
        sPData = add_CFO(pData,df,N);
        sSchData = add_CFO(schData,df,N);
        sSoumData = add_CFO(SoumTData,df,N);
        sLiuData = add_CFO(LiuData,df,N);
        sFangData = add_CFO(FangData,df,N);
        sRenData = add_CFO(RenData,df,N);
        sKimData = add_CFO(KimData,df,N);
        sZhangData = add_CFO(ZhangData,df,N);
        
        %Adding Additive White Gaussian Nois
        sPData = awgn(sPData,SNR(n),'measured');
        sSchData = awgn(sSchData,SNR(n),'measured');
        sSoumData = awgn(sSoumData,SNR(n),'measured');
        sLiuData = awgn(sLiuData,SNR(n),'measured');
        sRenData = awgn(sRenData,SNR(n),'measured');
        sFangData = awgn(sFangData,SNR(n),'measured');
        sKimData = awgn(sKimData,SNR(n),'measured');
        sZhangData = awgn(sZhangData,SNR(n),'measured');
        
        %Estimating and revising the fraction frequency offset(FFO)
        ffo = pAlgo.FFOEstimator(N,sPData(1,Ng+1:Ns));%Using the proposed FFO Algo to Estimate;
        %ffo = 0;       
        sPData = FFO_revise(sPData,ffo,N);
        sSchData = FFO_revise(sSchData,ffo,N);
        sSoumData = FFO_revise(sSoumData,ffo,N);
        sLiuData = FFO_revise(sLiuData,ffo,N);
        sRenData = FFO_revise(sRenData,ffo,N);
        sFangData = FFO_revise(sFangData,ffo,N);
        sKimData = FFO_revise(sKimData,ffo,N);
        sZhangData = FFO_revise(sZhangData,ffo,N);
        %Estimating Integer Frequency offset(IFO)
        pr = sPData(1,Ng+1:Ns); %remove the cyclic prefix(CP)
        pIFO = pAlgo.IFOEstimator(pr,N,pTS);
        if pIFO == ifo_df
            pSR(n) = pSR(n) + 1;
        end
        
        sSchData_r1 = sSchData(1,Ng+1:Ns) ;
        sSchData_r2 = sSchData(1,Ns+Ng+1:2*Ns);
        SchIFO = schAlgo.IFOEstimator(sSchData_r1,sSchData_r2,N,Sch_v);
        if SchIFO == ifo_df
            schSR(n) = schSR(n) + 1;
        end
        
        Soum_r = sSoumData(1,Ng+1:Ns);
        SoumIFO = SoumAlgo.IFOEstimator(Soum_r,N,SoumFData);
        if SoumIFO == ifo_df
           soumSR(n) = soumSR(n) + 1;
        end
        
        RenIFO = RAlgo.IFOEstimator(sRenData(1,Ng+1:Ns),N,RenTS);
        if RenIFO == ifo_df
            renSR(n) = renSR(n) + 1;
        end
        
        LiuIFO = LiuAlgo.IFOEstimator(sLiuData(1,Ng+1:Ns),N,LiuTS);
        if LiuIFO == ifo_df
            liuSR(n) = liuSR(n) + 1;
        end   
        
        FangIFO = FangAlgo.IFOEstimator(sFangData(1,Ng+1:Ns),N,Fang_v,Fang_mu);
        if FangIFO == ifo_df
            fangSR(n) = fangSR(n) + 1;
        end
 
        KimIFO = KimAlgo.IFOEstimator(sKimData(1,Ng+1:Ns),N,KimDiffSeq);
        if KimIFO == ifo_df
            kimSR(n) = kimSR(n) + 1;
        end    
        
        ZhangIFO = ZhangAlgo.IFOEstimator(sZhangData(1,Ng+1:Ns),N,ZhangDiffSeq);
        if ZhangIFO == ifo_df
            zhangSR(n) = zhangSR(n) + 1;
        end  
    end
    pSR(n) = pSR(n)/simu_times;
    schSR(n) = schSR(n)/simu_times;
    soumSR(n) = soumSR(n)/simu_times;
    fangSR(n) = fangSR(n)/simu_times;
    renSR(n) = renSR(n)/simu_times;
    liuSR(n) = liuSR(n)/simu_times;
    kimSR(n) = kimSR(n)/simu_times;
    zhangSR(n) = zhangSR(n)/simu_times;
end
figure 
d = 1:length(SNR);
semilogy(config.SNR,1-schSR(d),'-o','LineWidth',1);
hold on
semilogy(config.SNR,1-soumSR(d),'-db','LineWidth',1);
hold on
semilogy(config.SNR,1-renSR(d),'-x','LineWidth',1);
hold on
semilogy(config.SNR,1-liuSR(d),'-v','LineWidth',1);
hold on
semilogy(config.SNR,1-fangSR(d),'-s','LineWidth',1);
hold on
semilogy(config.SNR,1-kimSR(d),'-h','LineWidth',1);
hold on
semilogy(config.SNR,1-zhangSR(d),'-*','LineWidth',1);
hold on
semilogy(config.SNR,1-pSR(d),'-^k','LineWidth',1);
hold on
xlabel('SNR(dB)'); 
ylabel('Probabilty of Failure (POF)'); 
legend('Schmidl','Soumitra','Ren','Liu Chun','FangYibo','YHKim','Zhangjie','Proposed');
figure
d = 1:length(SNR);
plot(config.SNR,schSR(d),'-oc','LineWidth',1);
hold on
plot(config.SNR,soumSR(d),'-db','LineWidth',1);
hold on
plot(config.SNR,renSR(d),'-x','LineWidth',1);
hold on
plot(config.SNR,liuSR(d),'-v','LineWidth',1);
hold on
plot(config.SNR,fangSR(d),'-s','LineWidth',1);
hold on
plot(config.SNR,kimSR(d),'-h','LineWidth',1);
hold on
plot(config.SNR,zhangSR(d),'-*','LineWidth',1);
hold on
plot(config.SNR,pSR(d),'-^k','LineWidth',1);
hold on
xlabel('SNR(dB)'); 
ylabel('Probabilty of Success (POS)'); 
legend('Schmidl','Soumitra','Ren','Liu Chun','FangYibo','YHKim','Zhangjie','Proposed');

function[recv_sig] = add_CFO(trans_sig,df,N)
    %===========Ìí¼ÓÆµÆ«===========
    recv_sig = zeros(1,length(trans_sig));
    parfor k=1:length(trans_sig)
        %Ìí¼ÓÆµÆ«
        recv_sig(k) = trans_sig(k)*exp(sqrt(-1)*2*pi*df*k/N);
    end  
end

function [revised_recv] = FFO_revise(recv,df,N)
    parfor k=1:length(recv)
        recv(k) = recv(k)*exp(-sqrt(-1)*2*pi*df*k/N);
    end
    revised_recv = recv;    
end