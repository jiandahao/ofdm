classdef YHKimAlgo
    methods
        function [data, DiffSeq] = generate_data(obj,N,Ng)
        %============Construct Data (Same as Schimidl method)=================
            %seq = randi([0,15],1,N/2);
            QAMTable =[7+7i,-7+7i,-7-7i,7-7i]; 
            buf=QAMTable(randi([0,3],N/2,1)+1); 
           % buf = qammod(seq,16);
            Fdata=zeros(1,N); 
            index = 1; 
            for n=1:2:N 
                 Fdata(n)=sqrt(2)*buf(index); 
                 index=index+1; 
            end
            DiffSeq = zeros(1,N);
            for k = 1:2:N-2
                DiffSeq(k) = Fdata(k)/Fdata(k+2);%Fdata(2*k + 1)/Fdata(2*k+3);
            end
            data = ifft(Fdata,N);   %[A A] 
            data = [data(1,N-Ng+1:N) data]; 
        end
        
        function [df] = IFOEstimator(obj,recv_sig,N,DiffSeq)
            Y = fft(recv_sig);
            R = zeros(1,N/4+1);
            E = zeros(1,N/4+1);
            for g = 0:N/4
                for k = 0:N/2-2
                    m1 = mod(2*k+2*g,N)+1;
                    m2 = mod(2*k+2*g+2,N)+1;
                    R(g+1) = R(g+1) + Y(m1)*conj(DiffSeq(2*k+1))*conj(Y(m2));
                    E(g+1) = E(g+1) + abs(Y(m1))^2;
                end
            end
            
            F = (abs(R).^2) ./(E.^2);
            [~,df] = max(F);
            df = 2*(df-1);
        end
    end
end