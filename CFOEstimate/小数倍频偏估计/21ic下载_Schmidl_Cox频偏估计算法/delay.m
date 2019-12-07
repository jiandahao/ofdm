% Program 2-7
% delay.m
% Gives delay to input signal
% 
% Programmed by H.Harada	 	
% 

function [iout,qout] = delay( idata, qdata , nsamp , idel )

%****************** variables *************************
% idata  input Ich data     
% qdata  input Qch data     
% iout   output Ich data
% qout   output Qch data
% nsamp   Number of samples to be simulated 
% idel   Number of samples to be delayed
%******************************************************

iout=zeros(1,nsamp);
qout=zeros(1,nsamp);

if idel ~= 0 
  iout(1:idel) = zeros(1,idel);
  qout(1:idel) = zeros(1,idel);
end

iout(idel+1:nsamp) = idata(1:nsamp-idel);
qout(idel+1:nsamp) = qdata(1:nsamp-idel);

% ************************end of file***********************************
