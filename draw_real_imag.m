N=256;   
QAMTable=[7+7i,-7+7i,-7-7i,7-7i];
src = QAMTable(randi([0,3],N,1)+1); 
sym = ifft(src); 
figure(1);
subplot(1,2,1);
plot(real(sym));
xlabel('QAMTable:real after ifft'); 
subplot(1,2,2);
plot(imag(sym));
xlabel('QAMTable:imag after ifft'); 

a = zeros(1,N);
for n=0:N-1
    %a(n+1) = 7*sqrt(2)*exp(sqrt(-1)*pi*n*n/N);%orignal
    a(n+1) = exp(sqrt(-1)*pi*n*n/N);
end
signal = [sym a sym];
figure(2);
subplot(1,2,1);
plot(real(signal));
xlabel('CAZAC:real'); 
subplot(1,2,2);
plot(imag(signal));
xlabel('CAZAC:imag'); 

a = zeros(1,N);
for n=0:N-1
    %a(n+1) = 7*sqrt(2)*exp(sqrt(-1)*pi*n*n/N);%orignal
    a(n+1) = exp(sqrt(-1)*pi*n*n/N);
end
signal = [sym ifft(a) sym];
figure(3);
subplot(1,2,1);
plot(real(signal));
xlabel('CAZAC:real after ifft'); 
subplot(1,2,2);
plot(imag(signal));
xlabel('CAZAC:imag after ifft'); 

a= zeros(1,N);
for n=0:N-1
    a(n+1) = 7*sqrt(2)*exp(sqrt(-1)*pi*n*n/N);%orignal
    %a(n+1) = exp(sqrt(-1)*pi*n*n/N);
end
signal = [sym a sym];
figure(4);
subplot(1,2,1);
plot(real(signal));
xlabel('7*sqrt(2)*CAZAC:real '); 
subplot(1,2,2);
plot(imag(signal));
xlabel('7*sqrt(2)*CAZAC:imag'); 

a= zeros(1,N);
for n=0:N-1
    a(n+1) = 7*sqrt(2)*exp(sqrt(-1)*pi*n*n/N);%orignal
    %a(n+1) = exp(sqrt(-1)*pi*n*n/N);
end
b = ifft(a);
signal = [sym b sym];
figure(5);
subplot(1,2,1);
plot(real(signal));
xlabel('7*sqrt(2)*CAZAC:real after ifft'); 
subplot(1,2,2);
plot(imag(signal));
xlabel('7*sqrt(2)*CAZAC:imag after ifft'); 


