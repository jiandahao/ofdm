classdef ZhouLuAlgo
    methods
        function [data,TS] =  generate_data(obj,N,Ng)
           QAMTable =[7+7i,-7+7i,-7-7i,7-7i]; 
            buf=QAMTable(randi([0,3],N/2,1)+1); 
           % buf = qammod(seq,16);
            Fdata=zeros(1,N); 
            index = 1; 
            for n=1:2:N 
                 Fdata(n)=sqrt(2)*buf(index); 
                 index=index+1; 
            end
            TS = ifft(Fdata,N);   %[A A] 
            data = [TS(1,N-Ng+1:N) TS];             
        end
        
        function [df] = IFOEstimator(obj,r,N,TS)
            P = zeros(1,N/2+1);
            for g = 0:N/2
                for k=0:N-1
                    P(g+1) = P(g+1) + abs(r(k+1) - exp(2*pi*g/N)*TS(k+1));
                end
            end
            [~,df] = max(P);
            df = df-1;
        end
    end
end