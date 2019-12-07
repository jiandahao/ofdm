function [ ] = draw( name,M,n )
%UNTITLED9 此处显示有关此函数的摘要
%   此处显示详细说明
figure(n);
d=1:1:400;
% range=-177:1:223;
if strcmp(name,'Park')
   plot(d,M(d+128));
    h1 =legend([name,'''s method']); 
else if strcmp(name,'Proposed')
        plot(d,M(d));
        h1=legend([name,' method']);
    else
%[a b]=max(M);
    plot(d,M(d));
    h1=legend([name,'''s method']);
    end
end
h2=xlabel('Timing sampling'); 
h3=ylabel('Correlation Amplitude'); 
axis([0,400,0,1.1]);
set(gca,'XTick',0:25:400);
set(gca,'YTick',0:0.1:1.1);
set(h1,'Fontsize',18);
set(h2,'Fontsize',14);
set(h3,'Fontsize',14);
% if strcmp(name,'Park')
%    plot(range,M(d+128));
%     h =legend([name,'''s method']); 
% else if strcmp(name,'Proposed')
%         plot(range,M(d));
%         h=legend([name,' method']);
%     else
% %[a b]=max(M);
%     plot(range,M(d));
%     h=legend([name,'''s method']);
%     end
% end
% xlabel('Time (sample)'); 
% ylabel('Timing Metric'); 
% axis([-200,200,0,1.1]); 
% set(h,'Fontsize',18)
end

