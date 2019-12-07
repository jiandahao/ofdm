function [ trasmit_data,args ] = generate_trasmit_data( N,Ng,predata,suffixdata,type)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
%Ns = N + Ng;
QAMTable=[7+7i,-7+7i,-7-7i,7-7i]; 
if strcmp(type,'ren')
    a = zeros(1,N/2);
    for n=1:N/2
       % a(n) = 7*sqrt(2)*exp(2*sqrt(-1)*pi*n*n/N);%original
       a(n) =  7*sqrt(2)*exp(2*sqrt(-1)*pi*n*n/N);
    end
    A = ifft(a);
    %A=a;
    pn= 2*(rand(1,N)>0.5)-1;
    signal = pn.*[A A];
    cp_train = [signal(1,N-Ng+1:N) signal];
    trasmit_data = [predata cp_train suffixdata]; 
    args = pn;
end
if strcmp(type,'fang')
    a = zeros(1,N/2);
    for n=1:N/2
       % a(n) = 7*sqrt(2)*exp(2*sqrt(-1)*pi*n*n/N);%original
        a(n) =  7*sqrt(2)*exp(2*sqrt(-1)*pi*n*n/N);
    end
    A = ifft(a);
    %A=a;
    m = 1.2*(rand(1,N/2))-0.2;
    v= exp((sqrt(-1)*pi).*m);
    signal = [v.*A A];
    cp_train = [signal(1,N-Ng+1:N) signal];
    trasmit_data = [predata cp_train suffixdata];  
    args = v;
end
if strcmp(type,'shao')%按照原论文的方法，CAZAC是不经过IFFT处理的
    a = zeros(1,N/4);
    for n=0:N/4-1
        %a(n+1) = 7*sqrt(2)*exp(sqrt(-1)*pi*n*n/N);%orignal
        a(n+1) = exp(sqrt(-1)*pi*n*n/N);
    end
    %A = ifft(a);
    A=a;
    B = conj(a(1,N/4:-1:1));
    signal = [A B conj(A) conj(B)];
    cp_train = [signal(1,N-Ng+1:N) signal];
    trasmit_data = [predata cp_train suffixdata];
end
if strcmp(type,'wang')
    a = zeros(1,N/4);
    for n=1:N/4
        a(n) = 7*sqrt(2)*exp(4*sqrt(-1)*pi*n*n/N);
    end
    A = ifft(a);
    pn= 2*(rand(1,N/4)>0.5)-1;
    B = pn.*A;
    signal =[A B A B];
    cp_train = [signal(1,N-Ng+1:N) signal];
    trasmit_data = [predata cp_train suffixdata]; 
    args = pn;
end
if strcmp(type,'liubin')
    a = zeros(1,N/4);
    for n=1:N/4
        a(n) = 7*sqrt(2)*exp(4*sqrt(-1)*pi*n*n/N);%sqrt(2)*
    end
    % A = ifft(a);
    A=a;
    B = conj(a(1,N/4:-1:1));
    pn= 2*(rand(1,N/2)>0.5)-1;
    C = pn.*[A B];
    signal = [C A B];
    cp_train = [signal(1,N-Ng+1:N) signal];
    trasmit_data = [predata cp_train suffixdata];
    args = pn;
end
if strcmp(type,'park')
    pn=rand(1,N/4)>0.5;
    x=ifft(pn*7*sqrt(2));
    sch = [x x(1,N/4:-1:1) conj(x) conj(x(1,N/4:-1:1))];
    cp_train = [sch(1,N-Ng+1:N) sch];
    trasmit_data =[predata cp_train suffixdata];          
end

if strcmp(type,'cazac')
    a = zeros(1,N/4);
%     for n=1:N/4
%         a(n) = 7*sqrt(2)*exp(4*sqrt(-1)*pi*n*n/N);
%     end
    for n=0:N/4-1
       % a(n+1) =7*sqrt(2)*exp(4*sqrt(-1)*pi*n*n/N);%original
       a(n+1) =exp(4*sqrt(-1)*pi*n*n/N);
    end
    %A = ifft(a);
    A =a;
    B = conj(A(1,N/4:-1:1));
    C= zeros(1,N/4);
    for n=1:1:N/4
        if mod(n,2)
           % C(n) = (-1)*conj(B(n));
           C(n) = (-1)*B(n);
        else
            %C(n) = conj(B(n));
            C(n) = B(n);
        end
    end
  % D = conj(C); 
    D = C;
    signal = [A B C D];
    cp_train = [signal(1,N-Ng+1:N) signal];
    trasmit_data = [predata cp_train suffixdata];
end
if strcmp(type,'schmidl')
    buf=QAMTable(randi([0,3],N/2,1)+1); %1x(N/2)的矩阵
    x=zeros(1,N); 
    index = 1; 
    for n=1:2:N 
         x(n)=buf(index); 
         index=index+1; 
    end; 
    sch = ifft(x);   %[A A]的形式 
    cp_train = [sch(1,N-Ng+1:N) sch];
    trasmit_data = [ predata cp_train suffixdata]; 
end
if strcmp(type,'minn')
    buf=QAMTable(randi([0,3],N/2,1)+1);
    x=zeros(1,N/2 ); 
    index = 1; 
    for n=1:2:N/2 
         x(n)=buf(index); 
         index=index+1; 
    end; 
    sch = ifft(x);   %[A A]的形式 
    sch2=[sch (-1).*sch];%[A A -A -A]形式
    cp_train = [sch2(1,N-Ng+1:N) sch2];
    trasmit_data = [predata cp_train suffixdata];  
    end

end

