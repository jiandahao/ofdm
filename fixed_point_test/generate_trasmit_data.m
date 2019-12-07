function [ trasmit_data,args ] = generate_trasmit_data( N,Ng,predata,suffixdata)

    a = zeros(1,N/4);

    for n=0:N/4-1
       a(n+1) =exp(4*sqrt(-1)*pi*n*n/N);
    end
    A =a;
    B = conj(A(1,N/4:-1:1));
    C= zeros(1,N/4);
    for n=1:1:N/4
        if mod(n,2)
            C(n) = (-1)*conj(B(n));
        else
            C(n) = conj(B(n));
        end
    end
   D = conj(C); 
    signal = [A B C D];
    cp_train = [signal(1,N-Ng+1:N) signal];
    trasmit_data = [predata cp_train suffixdata];



end

