% @inproceedings{Liu2010A,
%   title={A robust frequency offset estimation scheme for OFDM systems},
%   author={Liu, Chun Guo and Li, Li Zhong and Zhang, Xiao Yong},
%   booktitle={Global Mobile Congress},
%   year={2010},
% }
classdef LiuChunAlgo
   methods 
       function [data,TS] = generate_data(obj,N,Ng)
        %============Construct Data=================
            QAMTable =[7+7i,-7+7i,-7-7i,7-7i]; 
            buf=QAMTable(randi([0,3],N/2,1)+1); 
            Fdata=zeros(1,N); 
            index = 1; 
            for n=1:2:N 
                 Fdata(n)=buf(index); 
                 index=index+1; 
            end; 
            data = ifft(Fdata,N);   %[A A] 
            TS = Fdata;
            data = [data(1,N-Ng+1:N) data];
       end
       
       function [df] = IFOEstimator(obj,recv_sig,N,TS)
           X = fft(recv_sig);
           R = zeros(1,N+1);
           for g = -N/2+1:N/2-1
               for k=0:N-1
                   R(g+N/2) = R(g+N/2) + X(mod(k+g,N)+1)*conj(TS(k+1)); 
               end
           end
           E1 = 0;
           E2 = 0;
           for k =0:N-1
               E1 = E1 + abs(X(k+1))^2;
               E2 = E2 + abs(TS(k+1))^2;
           end
           H = (abs(R).^2)./(E1*E2);
           [~,df] = max(H);
           df = df - N/2;
       end
   end
end