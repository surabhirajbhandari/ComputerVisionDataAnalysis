



clear; clc;
filename = "Rajbhandari S-Project-W2425.xlsx";
tableIn = xlsread(filename);
x= tableIn(1:70)
y=tableIn(71:130)

%Target Absent
figure(1)
[booth0stat,booth0sam] = bootstrp(5000,[],x);
y1 = x(booth0sam);
nb=5000;% number boot samples
ym=bootstrp(nb,@mean,x);%
histfit(ym,18,"Weibull") 
hold on
y25=prctile(ym,2.5); % 2.5%
y975=prctile(ym,97.5); % 97.5%

xlabel('boot samples value'),ylabel('frequency')
tit1={['Number of boot samples = ', num2str(nb)],['µ = ',num2str(mean(ym)),',σ = ',num2str(std(ym))],['95% CI = [',num2str(y25),', ',num2str(y975),']']};
dim = [0.13 0.62 0.3 0.3];
annotation('textbox',dim,'String',tit1,'FitBoxToText','on')
title("Statistics of Mean (H_0): Mean(µ), Std. Dev.(σ), 95% CI")

%Target Present
figure(2)
[booth1stat,booth1sam] = bootstrp(5000,[],y);
y2 = y(booth1sam);
yn=bootstrp(nb,@mean,y);%
histfit(yn,18,"Rician") 
hold on
y225=prctile(yn,2.5); % 2.5%
y2975=prctile(yn,97.5); % 97.5%

xlabel('boot samples value'),ylabel('frequency')
t={['Number of boot samples = ', num2str(nb)],['µ = ',num2str(mean(yn)),', σ = ',num2str(std(yn))],['95% conf. Int = [',num2str(y225),', ',num2str(y2975),']']};
annotation('textbox',dim,'String',t,'FitBoxToText','on')
title("Statistics of Mean (H_1): Mean(µ), Std. Dev.(σ), 95% CI")

% AUC
N0=length(x)
N1=length(y)
dat=[x;y]
resp = [zeros(N0,1); ones(N1,1)];

[pf,pd,T,AUC]=perfcurve(resp,dat,1);
btarea=bootstrp(nb,@calculateROCarea,resp,dat);
figure(3), histfit(btarea,20)
y325=prctile(btarea,2.5); 
y3975=prctile(btarea,97.5); 
xlabel('boot samples value'),ylabel('frequency')
tit1={['Number of boot samples = ', num2str(nb)],['µ = ',num2str(mean(btarea)),', σ = ',num2str(std(btarea))],['95% CI = [',num2str(y325),', ',num2str(y3975),']']};
annotation('textbox',dim,'String',tit1,'FitBoxToText','on')
title("Statistics of AUC: Mean(µ), Std. Dev.(σ), 95% CI")

function area=calculateROCarea(resp,dat)
[pf,pd,T,AUC]=perfcurve(resp,dat,1);
area=AUC;
end
