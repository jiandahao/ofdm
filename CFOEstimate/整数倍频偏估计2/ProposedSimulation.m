clear all;
clc

config = OFDMSystemConfig;
N = config.N;
Ng = config.Ng;

%DrawMetricUnderDiffFFO(N,Ng);
%DrawMetricUnderDiffIFO(N,Ng);
%DrawMetricUnderDiffSNR(N,Ng);
CalculateProperW(N,Ng);

function DrawMetricUnderDiffSNR(N,Ng)
    W = N/8;
    snr = [-10:2:10];
    rangeSize = 2*W+1;
    M = zeros(length(snr),rangeSize);
    estimateVal = zeros(1,length(snr));
    [pData,TS] =  generate_data(N,Ng);
    %ifo = 20;
    Ns=N+Ng;
    for n = 1:length(snr)
        spData = add_CFO(pData,16,N);
        spData = awgn(spData,snr(n),'measured');
        [M(n,1:rangeSize),estimateVal(n)] = get_IFO(spData(1,Ng+1:Ns),N,TS);
    end
    figure
    for n = 1:length(snr)
        x = estimateVal(n)-W:estimateVal(n)+W;
        y = snr(n)*ones(1,rangeSize);
        plot3(y,x,M(n,1:rangeSize),'black');
        hold on;
    end
    x1 = xlabel('SNR(dB)'); 
    y1 = ylabel('IFO Scan Window');
%      set(gca,'xtick',[-8:4:8]);
%      set(x1,'Rotation',15);
%      set(y1,'Rotation',-25);
%    set(gca,'ytick',[-16:4:16]);
    zlabel('Metrix curve of the Correlation');
    grid on
end

function DrawMetricUnderDiffIFO(N,Ng)
    W = N/8;
    ifo = [-8:2:8];
    rangeSize = 2*W+1;
    M = zeros(length(ifo),rangeSize);
    estimateVal = zeros(1,length(ifo));
    [pData,TS] =  generate_data(N,Ng);
    %ifo = 20;
    Ns=N+Ng;
    for n = 1:length(ifo)
        spData = add_CFO(pData,ifo(n),N);
        spData = awgn(spData,15,'measured');
        [M(n,1:rangeSize),estimateVal(n)] = get_IFO(spData(1,Ng+1:Ns),N,TS);
    end
    figure
   
    for n = 1:length(ifo)
        x = estimateVal(n)-W:estimateVal(n)+W;
        y = ifo(n)*ones(1,rangeSize);
        plot3(y,x,M(n,1:rangeSize),'black');
        hold on;
    end
    x1 = xlabel('IFO'); 
    y1 = ylabel('IFO Scan Window');
     set(gca,'xtick',[-8:4:8]);
     set(x1,'Rotation',15);
     set(y1,'Rotation',-25);
%    set(gca,'ytick',[-16:4:16]);
    zlabel('Metrix curve of the Correlation');
    grid on
end

function DrawMetricUnderDiffFFO(N,Ng)
    W = N/2;
    ffo = -0.5:0.1:0.5;
    rangeSize = 2*W+1;
    M = zeros(length(ffo),rangeSize);
    estimateVal = zeros(1,length(ffo));
    [pData,TS] =  generate_data(N,Ng);
    ifo = 20;
    Ns=N+Ng;
    for n = 1:length(ffo)
        spData = add_CFO(pData,ifo+ffo(n),N);
        spData = awgn(spData,15,'measured');
        [M(n,1:rangeSize),estimateVal(n)] = get_IFO(spData(1,Ng+1:Ns),N,TS);
    end
    figure
    for n = 1:length(ffo)
        x = estimateVal(n)-W:estimateVal(n)+W;
        y = ffo(n)*ones(1,rangeSize);
        plot3(y,x,M(n,1:rangeSize),'black');
        hold on;
    end
    x1 = xlabel('FFO'); 
    y1 = ylabel('IFO');
    zlabel('Metrix curve of the Correlation');
    grid on
end

function [H,val] = get_IFO(recv_sig,N,TS)
    Q=0;
    %==============´Ö¹À¼Æ=============
    for k=1:N-1
       Q1 = recv_sig(k)*conj(TS(k));
       Q2 = recv_sig(k+1)*conj(TS(k+1));
       Q = Q + conj(Q1)*Q2;
    end
    idf_coarse = N*angle(Q)/(2*pi);
    if idf_coarse>=0
        idf = floor(idf_coarse);
    else
        idf = ceil(idf_coarse);
    end
    val = idf;
    %==============Ï¸¹À¼Æ==============
    W= N/8;
    P = zeros(1,2*W+1);
    E = 0;
    for g = idf-W:idf+W
        for n=1:N
            P(g-idf+W+1) = P(g-idf+W+1) + recv_sig(n)*conj(TS(n))*exp(-sqrt(-1)*2*pi*g*n/N);
        end
    end

    for n=1:N
        E = E+abs(recv_sig(n))^2;
    end

     H = abs(P).^2/E;
%     [~,df]= max(H);
%     df = df +idf-W-1;
end

function [data,TS] = generate_data(N,Ng)
    a = zeros(1,N/4);
    mu = N/4-1;
    for n=0:N/4-1
       a(n+1) =exp(4*sqrt(-1)*mu*pi*n*n/N);
    end
    A = ifft(a);
    B = conj(A(1,N/4:-1:1));
    C= zeros(1,N/4);
    for n=1:1:N/4 
        if mod(n,2)
           C(n) = (-1)*B(n);
        else
           C(n) = B(n);
        end
    end
    D=C;
    data = [A B C D];
    TS = data;
    data =[data(1,N-Ng+1:N) data];%Ê±ÓòÉÏµÄÑµÁ··ûºÅ
end

function CalculateProperW(N,Ng)
    ifo = 16;
    W = N/8; % set a initial value;
    snr= [-10:2:10];
    simu_times = 1000;
    Ns = N + Ng;
    [pData,TS] =  generate_data(N,Ng);
    maxDw = zeros(1,length(snr)); % store the max offset between Estimated ifo and real ifo;
    avergeDw = zeros(1,length(snr));
    parfor n = 1:length(snr)
        for s = 1:simu_times
        ifo = randi([0,32]);
        ffo = -0.5+rand();
        df = ifo + ffo;
        spData = add_CFO(pData,df,N);
        spData = awgn(spData,snr(n),'measured');
        ffo = FFOAlgoProposed(N,spData);
        spData = FFO_revise(spData,ffo,N);
        
        [~,EstIFO] = get_IFO(spData(1,Ng+1:Ns),N,TS);
        
        maxDw(n) = max(maxDw(n),abs(EstIFO - ifo));
        avergeDw(n) =  avergeDw(n) + abs(EstIFO - ifo);
        end
        avergeDw(n) = avergeDw(n)/simu_times;
    end
    figure
    d = 1:length(snr);
    plot(snr,maxDw(d),'-*');
    hold on
    plot(snr,avergeDw(d),'-h');
    xlabel('SNR(dB)'); 
    ylabel('Estimated Integer Frequency Offset'); 
    legend('×î´ó¾ø¶ÔÆ«²î','Æ½¾ù¾ø¶ÔÆ«²î');
    grid on

end

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