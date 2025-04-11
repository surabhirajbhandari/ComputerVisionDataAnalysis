close all;
clear;
 
[status,sheets] = xlsfinfo('subu.xlsx'); 
[A,names,raw] =xlsread('subu.xlsx',1); 
 
A; 
Abs=A(1:70);
N0=length(Abs);
Prs=A(71:end); 
N1=length(Prs);
dat=[Abs;Prs];
resp =[zeros(N0,1);ones(N1,1)];
[pf,pd,T,AUC,opt]=perfcurve(resp,dat,1);
parg0=fitdist(Abs,'gamma')
parg1=fitdist(Prs,'gamma')
thr = 0:0.01:100;
aprs=parg1.a
aabs=parg0.a
bprs=parg1.b
babs=parg0.b
PF=1-cdf('gamma',thr,aabs,babs)
PD=1-cdf('gamma',thr,aprs,bprs)
AUCG=0.5 + polyarea(PF,PD)

xx=0:0.001:1.2*max(dat);

%{
figure(5)
histogram(Abs, 'normalization', 'pdf')
hold on
histogram(Prs, 'normalization', 'pdf')

plot (xx, fx1,'-r', 'linewidth',1.5)
plot (xx, fx2,'--g', 'linewidth',1.5)
intersect=find(round(fx1,2)==round(fx2,2))
finds max & median
peaks_fx1 = findpeaks(fx1);
peaks_fy1 = findpeaks(fx2);
plot(1.063, peaks_fx1, "red", "Marker","o", "LineWidth", 2)
plot(2.56, peaks_fy1, '*k', "Marker", "square", "LineWidth", 2)
%}
fx1=ksdensity(Abs,xx);
fx2=ksdensity(Prs,xx);

%Use Midpoint as Threshold
vt=1.8048;
CorrectTargetDet = Prs > vt;
NcM = sum(CorrectTargetDet);
Miss = Abs > vt;
NfM = sum(Miss);
N=N0+N1
%Confusion Matrix
confusionMM = [(N0 - NfM), (NfM); (N1 -NcM) (NcM)]
errorrateM = (NfM + (N1 - NcM))/N
ppvM = NcM/(NfM + NcM)
pmM = NfM/130;
pfM = NcM/130;
tranisitionmM = [(1 - pfM), pmM; pfM, (1-pmM)]



vt1 = 1.3987;
CorrectTargetDetI = Prs > vt1;
NcI = sum(CorrectTargetDetI);
MissI = Abs > vt1;
NfI = sum(MissI);
%Confusion Matrix
confusionMI = [(N0 - NfI), (NfI); (N1 -NcI) (NcI)]
errorrateI = (NfI + (N1 - NcI))/N
ppvI = NcI/(NfI + NcI)
pmI = NfI/130;
pfI = NcI/130;
tranisitionmI = [(1 - pfI), pmI; pfI, (1-pmI)]


%% data table
figure(1)
xlim([0,10])
ylim([0,5])
axis off
B=readmatrix('subu.xlsx','Sheet',1);
Abs_Sorted=B(1:70);
Prs_Sorted=B(71:end); 
meanAbs=mean(Abs_Sorted)
meanPrs=mean(Prs_Sorted)
varAbs=var(Abs_Sorted)
varPrs=var(Prs_Sorted)
pindex=(abs(meanAbs-meanPrs))/sqrt((varAbs+varPrs))
dats= reshape(Abs_Sorted,[10,7]);
text(1.5,2.89,'Target Absent (sorted)')
text(0.3,2.2,num2str(dats))
text(1,1.50,['mean µ_1 = ',num2str(meanAbs),','],'Color','b','FontWeight','bold')
text(2,1.50,['var σ_1^2 = ',num2str(varAbs)],'Color','b','FontWeight','bold')

ann=annotation('line',[.395 .395],[.34 .61])
ann.Color='r'
ann.LineWidth=1.5;

dats= reshape(Prs_Sorted,[10,6]);
text(4.1,2.89,'Target Present (sorted)')
text(4.3,2.2,num2str(dats))
text(4,1.50,['mean µ_2 = ',num2str(meanPrs)],'Color','b','FontWeight','bold')
text(5,1.50,['var σ_2^2 = ',num2str(varPrs)],'Color','b','FontWeight','bold')

title("Dervishi-Project_F2223")
set(get(gca,'title'),'Position',[3.4 3.3 1.00011])
text(2.5,0.4,['Performance Index PI = ',num2str(pindex)],'Color','b')

%% Confusion and Transition Matrices
disp(['Threshold (vt) = ',num2str(vt),' (midpoint)'])
disp(table(confusionMM,tranisitionmM,'VariableNames',{'Cx','Tx'}))
disp(['error rate = ', num2str(errorrateM),', ','PPV = ', num2str(ppvM)])

fprintf('\n')

disp(['Threshold (vt) = ',num2str(vt1),' (intersection)'])
disp(table(confusionMI,tranisitionmI,'VariableNames',{'Cx','Tx'}))
disp(['error rate = ', num2str(errorrateI),', ','PPV = ', num2str(ppvI)])

%% figure 2

dat=[Abs;Prs];
unsorted = [zeros(size(Abs));ones(size(Prs))];
unsorted = [unsorted, dat];
sorted = sortrows(unsorted, 2, 'descend');
sorted = [sorted; 0,0];
T = sorted(:, 2);
counts = zeros(size(sorted));
for i = 1:size(counts)
    Nc = 0;
    Nf = 0;
    if i == 0
        continue;
    else
        for j = 1:(i-1)
            if sorted(j, 1) == 1
                Nc = Nc + 1;
            else
                Nf = Nf + 1;
            end
        end
        counts(i, 1) = Nc;
        counts(i, 2) = Nf;
    end
end
P = zeros(size(counts));
for i = 1:111
    P(i, 1) = counts(i, 1)/N1;
    P(i, 2) = counts(i, 2)/N0;
end
dist = zeros(111,1);
for i = 1:111
    dist(i, 1) = pdist([0,1;P(i, 2), P(i, 1)],'euclidean');
end

[opt_dist, I] = min(dist);

opt_thresh = T(I);
%Confusion Matrix
CorrectTargetDeto = Prs > opt_thresh;
Nco = sum(CorrectTargetDeto);
Misso = Abs > opt_thresh;
Nfo = sum(Misso);
confusiono = [(N0 - Nfo), (Nfo); (N1 -Nco) (Nco)]
errorrateo = (Nfo + (N1 - Nco))/N
ppvo = Nco/(Nfo + Nco)
pmo = Nfo/130;
pfo = Nco/130;
tranisitionmo = [(1 - pfo), pmo; pfo, (1-pmo)]

figure(2);
xlim([0,1])
ylim([0,1])
plot(pf, pd,'k','LineWidth',1.5)
hold on

plot(PF,PD,':m','LineWidth',2)
plot(opt(1),opt(2), ".r", markersize=15)
plot([0,1],[0,1], '--g','LineWidth',1.5)
plot([0,opt(1)],[1,opt(2)],'-b','linewidth',1.25)
text(0.15,0.85,['d_o_p_t = ',num2str(opt_dist)])
text(opt(1)*1.1,opt(2)*.95,['[P_F = ', num2str(opt(1)),', P_D = ', num2str(opt(2)),']'])
xlabel('P_F = [1-Specificity]'),ylabel('Sensitivity = [1-P_M]')
strAUC=['Empirical ROC: AUC = ',num2str(AUC)]
strAUCG=['Bigamma fit: AUC = ',num2str(AUCG)]

legend(strAUC,strAUCG,'OOP (Neyman Pearson)','Location', 'Southeast')
title(['Optimal Threshold Value= ',num2str(opt_thresh)],'Color','k')


%% figure 3

figure(3)

xx=0:0.0001:1.2*max(Abs);
xx1 = xx;fx1=ksdensity(Abs,xx);
plot(xx,fx1,'k' ,'linewidth',1.5)
xlim([0,5])
xlabel('values(v)'),ylabel('estimated pdf')
hold on 
xx=0:0.0001:1.2*max(Prs);
fx2=ksdensity(Prs,xx); 
plot(xx,fx2,'--r','linewidth',1.5)
title('density fit')
hold on
text(.1,1.5,'PM')
idx_1 = find(abs(xx1 - opt_thresh) < 0.00001);
X1=[xx1(idx_1:end),fliplr(xx1(idx_1:end))];
Y1=[fx1(idx_1:end),fliplr(fx1(idx_1:end) * 0)]; 
h1 = fill(X1,Y1,'g');
set(h1,'facealpha',.4,'facecolor','#808080')
idx_2 = find(abs(xx - opt_thresh) < 0.00001);
X=[xx(1:idx_2),fliplr(xx(1:idx_2))];
Y=[fx2(1:idx_2),fliplr(fx2(1:idx_2) * 0)];
h2 = fill(X,Y,'b');
set(h2,'facealpha',.4,'facecolor','blue')
legend('f_v(v|H_0)', 'f_v(v|H_0)')
ann=annotation('arrow',[opt_thresh*.255 opt_thresh*.255], [0.11 0.85])
ann.Color="black"
text(opt_thresh+.5,0.05,'P_F')
text(opt_thresh-.5,0.05,'P_M')
ax=opt_thresh+.25
text(ax,max(fx1),['v_T (optimal) = ',num2str(opt_thresh)])
