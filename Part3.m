clear;
 
[status,sheets] = xlsfinfo('subi_8.xlsx'); 
[A,names,raw] =xlsread('subi_8.xlsx',1); 
 
A; 
Abs=A(1:70);
N0=length(Abs);
Prs=A(71:end); 
N1=length(Prs);
da='Gamma'
d='Rician'

parAbs=fitdist(Abs,da)
parPrs=fitdist(Prs,d)
RNDAbs = random(parAbs,[numel(Abs) 2])
RNDPrs = random(parPrs,[numel(Prs) 2]);
 
pts = linspace(0, 15, 100)';



%% Figure 2
figure(2)
tiledlayout(2,2)

%Original Data

B=readmatrix('subi_8.xlsx','Sheet',1);
Abs_Sort=B(1:70);
Prs_Sort=B(71:end);
Abs_Sorted=sort(RNDAbs(:,1));
Prs_Sorted=sort(RNDPrs(:,1))
meanAbs=mean(Abs_Sort)
meanPrs=mean(Prs_Sort)
varAbs=var(Abs_Sort)
varPrs=var(Prs_Sort)
PI=(abs(meanAbs-meanPrs))/sqrt((varAbs+varPrs))


fx = ksdensity(reshape(Abs, [numel(Abs) 1]), pts);
fy = ksdensity(reshape(Prs, [numel(Prs) 1]), pts);

nexttile
xlim([0,10]),ylim([0 0.8])
hold on
plot(pts,fx,'r' ,'linewidth',1.5)
xlabel('Data'),ylabel('Estimated PDF')
hold on 
plot(pts,fy,'--k' ,'linewidth',1.5)
hold on


title(["Original data",['Performance Index = ',num2str(PI)]]) 
legend('Target Absent', 'Target Present')


%Arithmetic Mean
am0=(RNDAbs(:, 1) + RNDAbs(:, 2))/2;
am1=(RNDPrs(:, 1) + RNDPrs(:, 2))/2;
PIA = (abs(mean(am0)-mean(am1)))/(sqrt(var(am0)+var(am1)))
fxA=ksdensity(am0,pts)
fyA=ksdensity(am1,pts)

nexttile
xlim([0,10]),ylim([0 1])
hold on
plot(pts,fxA,'r' ,'linewidth',1.5)
xlabel('Data'),ylabel('Estimated PDF')
hold on 
plot(pts,fyA,'--k' ,'linewidth',1.5)
title(["Arithmetic Mean",['Performance Index = ',num2str(PIA)]]) 
legend('Target Absent', 'Target Present')


%Maximum
max0=max(RNDAbs(:, 1), RNDAbs(:, 2));
max1=max(RNDPrs(:, 1), RNDPrs(:, 2));
PIM = (abs(mean(max0)-mean(max1)))/(sqrt(var(max0)+var(max1)))
fxM=ksdensity(max0,pts)
fyM=ksdensity(max1,pts)

nexttile
xlim([0,10]),ylim([0 0.8])
hold on
plot(pts,fxM,'r' ,'linewidth',1.5)
xlabel('Data'),ylabel('Estimated PDF')
hold on 
plot(pts,fyM,'--k' ,'linewidth',1.5)
title(["Maximum",['Performance Index = ',num2str(PIM)]]) 
legend('Target Absent', 'Target Present')


%Geometric Mean
gm0=sqrt(RNDAbs(:, 1).*RNDAbs(:, 2));
gm1=sqrt(RNDPrs(:, 1).*RNDPrs(:, 2));
PIG = (abs(mean(gm0)-mean(gm1)))/(sqrt(var(gm0)+var(gm1)))
fxG=ksdensity(gm0,pts)
fyG=ksdensity(gm1,pts)

nexttile
xlim([0,10]),ylim([0 0.8])
hold on
plot(pts,fxG,'r' ,'linewidth',1.5)
xlabel('Data'),ylabel('Estimated PDF')
hold on
plot(pts,fyG,'--k' ,'linewidth',1.5)
title(["Geometric Mean",['Performance Index = ',num2str(PIG)]]) 
legend('Target Absent', 'Target Present')

%% Figure 3
figure(3)
resp =[zeros(N0,1);ones(N1,1)];
dat = [RNDAbs(:,1);RNDPrs(:,1)]
dat1 = [Abs(:,1);Prs(:,1)]
[pf,pd,t,AUC]=perfcurve(resp,dat1,1);
[pfA,pdA,t, AUCA]=perfcurve(resp,[am0;am1],1);
[pfG,pdG,t,AUCG]=perfcurve(resp,[gm0;gm1],1);
[pfM,pdM,t,AUCM]=perfcurve(resp,[max0;max1],1);

xlim([0,1])
ylim([0,1])
plot(pf, pd,'magenta','LineWidth',1.5)
hold on
plot(pfA, pdA,'--r','LineWidth',1.5)
plot(pfG, pdG,'-.b','LineWidth',1.5)
plot(pfM, pdM,':k','LineWidth',1.5)

plot([0,1],[0,1], '--g','LineWidth',1.5)
xlabel('Probability of False Alarm'),ylabel('Probability Detection')

legend(['Original data: AUC = ',num2str(AUC)],['Arithmetic Mean of two: AUC = ',num2str(round(AUCA,3))],['Geometric Mean of two: AUC = ',num2str(round(AUCG,3))],['Maximum of two: AUC = ',num2str(round(AUCM,3))],'location','southeast')
title('Empirical ROC curves: Pre- and Post-processing','Color','k')



%% Figure 4
figure(4)
axis off, grid off

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
for i = 1:131
    P(i, 1) = counts(i, 1)/N1;
    P(i, 2) = counts(i, 2)/N0;
end
dist = zeros(111,1);
for i = 1:131
    dist(i, 1) = pdist([0,1;P(i, 2), P(i, 1)],'euclidean');
end

[opt_dist, I] = min(dist);

ot = T(I);

NcO = sum(Prs(:,1) > ot);
NfO = sum(Abs(:,1) > ot);

Nct = sum(Prs(:,1) > ot);
Nft = sum(Abs(:,1) > ot);

dat=[am0;am1];
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
for i = 1:131
    P(i, 1) = counts(i, 1)/N1;
    P(i, 2) = counts(i, 2)/N0;
end
dist = zeros(111,1);
for i = 1:131
    dist(i, 1) = pdist([0,1;P(i, 2), P(i, 1)],'euclidean');
end

[opt_distA, I] = min(dist);

otA = T(I);

NcA = sum(am1 > otA);
NfA = sum(am0 > otA);

dat=[max0;max1];
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
for i = 1:131
    P(i, 1) = counts(i, 1)/N1;
    P(i, 2) = counts(i, 2)/N0;
end
dist = zeros(111,1);
for i = 1:131
    dist(i, 1) = pdist([0,1;P(i, 2), P(i, 1)],'euclidean');
end

[opt_distM, I] = min(dist);

otM = T(I);

NcM = sum(max1 > otM);
NfM = sum(max0 > otM);

dat=[gm0;gm1];
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
for i = 1:131
    P(i, 1) = counts(i, 1)/N1;
    P(i, 2) = counts(i, 2)/N0;
end
dist = zeros(131,1);
for i = 1:111
    dist(i, 1) = pdist([0,1;P(i, 2), P(i, 1)],'euclidean');
end

[opt_distG, I] = min(dist);
otg = T(I);

NcG = sum(gm1 > otg);
NfG = sum(gm0 > otg);

ERRORRATEA = (NfA + (N1 - NcA))
PPVA = NcA/(NfA + NcA)

ERRORRATEM = (NfM + (N1 - NcM))
PPVM = NcM/(NfM + NcM)

ERRORRATEG = (NfG + (N1 - NcG))
PPVG = NcG/(NfG + NcG)

ERRORRATE = (Nft + (N1 - Nct))
PPV = Nct/(Nft + Nct)

t = title("Performance Improvement Metrics: Dual Diversity")
t.Color='b'
text(0.25,0.9,["Target", "Not Detected"])
text(0.45,0.9,["Target", "Detected"])
text(0.6,0.9,["Signal Processing", "Algorithm"])
text(0,0.78,'Target Absent')
text(0,0.72,'Target Present')
text(0,0.58,'Target Absent')
text(0,0.52,'Target Present')
text(0,0.38,'Target Absent')
text(0,0.32,'Target Present')
text(0,0.18,'Target Absent')
text(0,0.12,'Target Present')
text(0.6,0.80,['No Processing(input): PI = ',num2str(round(PI,4))],'Color','r')
text(0.6,0.75,['*Error rate = ', num2str(round(ERRORRATE,3)), '/130,   AUC = ',num2str(round(AUC,3))])
text(0.6,0.70,['*Positive Predictive Value = ',num2str(PPV)])
text(0.6,0.60,['Arithmetic Mean: PI = ',num2str(round(PIA,3))],'Color','r')
text(0.6,0.55,['*Error rate = ', num2str(round(ERRORRATEA,3)), '/130,   AUC = ',num2str(round(AUCA,3))])
text(0.6,0.50,['*Positive Predictive Value = ',num2str(round(PPVA,3))])
text(0.6,0.40,['Maximum: PI = ',num2str(round(PIM,3))],'Color','r')
text(0.6,0.35,['*Error rate = ', num2str(round(ERRORRATEM,3)), '/130,   AUC = ',num2str(round(AUCM,3))])
text(0.6,0.30,['*Positive Predictive Value = ',num2str(round(PPVM,3))])
text(0.6,0.20,['Geometric Mean: PI = ',num2str(round(PIG,3))],'Color','r')
text(0.6,0.15,['*Error rate = ', num2str(round(ERRORRATEG,3)), '/130,   AUC = ',num2str(round(AUCG,3))])
text(0.6,0.10,['*Positive Predictive Value = ',num2str(round(PPVG,3))])
text(0.3,0.78,num2str((N0-NfO)),'FontWeight','bold','Color','r')
text(0.45,0.78,num2str((NfO)),'FontWeight','bold','Color','r')
text(0.3,0.72,num2str((N1-NcO)),'FontWeight','bold','Color','r')
text(0.45,0.72,num2str((NcO)),'FontWeight','bold','Color','r')
text(0.3,0.58,num2str((N0-NfA)),'FontWeight','bold','Color','r')
text(0.45,0.58,num2str((NfA)),'FontWeight','bold','Color','r')
text(0.3,0.52,num2str((N1-NcA)),'FontWeight','bold','Color','r')
text(0.45,0.52,num2str((NcA)),'FontWeight','bold','Color','r')
text(0.3,0.38,num2str((N0-NfM)),'FontWeight','bold','Color','r')
text(0.45,0.38,num2str((NfM)),'FontWeight','bold','Color','r')
text(0.3,0.32,num2str((N1-NcM)),'FontWeight','bold','Color','r')
text(0.45,0.32,num2str((NcM)),'FontWeight','bold','Color','r')
text(0.3,0.18,num2str((N0-NfG)),'FontWeight','bold','Color','r')
text(0.45,0.18,num2str((NfG)),'FontWeight','bold','Color','r')
text(0.3,0.12,num2str((N1-NcG)),'FontWeight','bold','Color','r')
text(0.45,0.12,num2str((NcG)),'FontWeight','bold','Color','r')
text(-0.05,0.35,'Single Iteration','FontWeight','bold','Color','r','Rotation',90,'FontSize',11)

annotation('line',[0.32,0.95],[0.79,0.79],'linewidth',1.5)
annotation('line',[0.32,0.95],[0.63,0.63],'linewidth',1.5)
annotation('line',[0.32,0.95],[0.47,0.47],'linewidth',1.5)
annotation('line',[0.32,0.95],[0.31,0.31],'linewidth',1.5)
annotation('line',[0.32,0.95],[0.15,0.15],'linewidth',1.5)

text(0.25,0,["*Error rate and positive predictive value estimated"," based on the Neyman-Pearson Threshold"])

%% Figure 5

AMAUC = zeros(100,1);
GMAUC = zeros(100,1);
MXAUC = zeros(100,1);
AME = zeros(100,1);
GME = zeros(100,1);
MXE = zeros(100,1);
AMPPV = zeros(100,1);
GMPPV = zeros(100,1);
MXPPV = zeros(100,1);
AMPI = zeros(1001,1);
GMPI = zeros(100,1);
MXPI = zeros(100,1);
for x=1:100
    parAbs=fitdist(Abs,da)
    parPrs=fitdist(Prs,d)
    RNDAbs = random(parAbs,[numel(Abs) 2])
    RNDPrs = random(parPrs,[numel(Prs) 2]);
    
    am0=(RNDAbs(:, 1) + RNDAbs(:, 2))/2;
    am1=(RNDPrs(:, 1) + RNDPrs(:, 2))/2;
    max0=max(RNDAbs(:, 1), RNDAbs(:, 2));
    max1=max(RNDPrs(:, 1), RNDPrs(:, 2));
    gm0=sqrt(RNDAbs(:, 1).*RNDAbs(:, 2));
    gm1=sqrt(RNDPrs(:, 1).*RNDPrs(:, 2));
    
    [pfA,pdA,t, AUCA]=perfcurve(resp,[am0;am1],1);
    [pfG,pdG,t,AUCG]=perfcurve(resp,[gm0;gm1],1);
    [pfM,pdM,t,AUCM]=perfcurve(resp,[max0;max1],1);
    
    dat=[am0;am1];
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
    for i = 1:131
        P(i, 1) = counts(i, 1)/N1;
        P(i, 2) = counts(i, 2)/N0;
    end
    dist = zeros(111,1);
    for i = 1:131
        dist(i, 1) = pdist([0,1;P(i, 2), P(i, 1)],'euclidean');
    end
    [opt_distA, I] = min(dist);
    otA = T(I);
    
    
    dat=[max0;max1];
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
    for i = 1:131
        P(i, 1) = counts(i, 1)/N1;
        P(i, 2) = counts(i, 2)/N0;
    end
    dist = zeros(111,1);
    for i = 1:131
        dist(i, 1) = pdist([0,1;P(i, 2), P(i, 1)],'euclidean');
    end
    [opt_distM, I] = min(dist);
    otM = T(I);
    
    
    dat=[gm0;gm1];
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
    for i = 1:131
        P(i, 1) = counts(i, 1)/N1;
        P(i, 2) = counts(i, 2)/N0;
    end
    dist = zeros(111,1);
    for i = 1:131
        dist(i, 1) = pdist([0,1;P(i, 2), P(i, 1)],'euclidean');
    end
    [opt_distG, I] = min(dist);
    otG = T(I);
    
    NcA = sum(am1 > otA);
    NfA = sum(am0 > otA);
    NcM = sum(max1 > otM);
    NfM = sum(max0 > otM);
    NcG = sum(gm1 > otG);
    NfG = sum(gm0 > otG);
    
    ERRORRATEA = (NfA + (N1 - NcA))
    PPVA = NcA/(NfA + NcA)
    ERRORRATEM = (NfM + (N1 - NcM))
    PPVM = NcM/(NfM + NcM)
    ERRORRATEG = (NfG + (N1 - NcG))
    PPVG = NcG/(NfG + NcG)
    
    PIA = (abs(mean(am0)-mean(am1)))/(sqrt(var(am0)+var(am1)))
    PIG = (abs(mean(gm0)-mean(gm1)))/(sqrt(var(gm0)+var(gm1)))
    PIM = (abs(mean(max0)-mean(max1)))/(sqrt(var(max0)+var(max1)))
    
    AMAUC(x)=AUCA
    MXAUC(x)=AUCM
    GMAUC(x)=AUCG
    AME(x)=ERRORRATEA
    GME(x)=ERRORRATEG
    MXE(x)=ERRORRATEM
    AMPPV(x)=PPVA
    GMPPV(x)=PPVG
    MXPPV(x)=PPVM
    AMPI(x)=PIA
    GMPI(x)=PIG
    MXPI(x)=PIM
    
end
%%
figure(5)
axis off, grid off
t = title(["Performance Improvement from Signal Processing(Dual Diversity)","Mean and Std. Deviation of Metrics: 100 Iterations"])
t.Color='b'

text(0.13,0.9,"Orig Data",'FontWeight','bold')
text(0.35,0.9,"AM",'FontWeight','bold')
text(0.55,0.9,"MX",'FontWeight','bold')
text(0.75,0.9,"GM",'FontWeight','bold')
text(-0.099,0.8,'AUC','Color','r','FontWeight','bold')
text(-0.099,0.6,["Error Counts","[out of 130]"],'Color','r','FontWeight','bold')
text(-0.099,0.4,'PPV ','Color','r','FontWeight','bold')
text(-0.099,0.2,'Perf. Index ','Color','r','FontWeight','bold')
text(0.9,0.84,'(mean)','Color','r','FontWeight','bold')
text(0.9,0.77,'(std. dev.)','Color','r','FontWeight','bold')
text(0.9,0.64,'(mean)','Color','r','FontWeight','bold')
text(0.9,0.57,'(std. dev.)','Color','r','FontWeight','bold')
text(0.9,0.44,'(mean)','Color','r','FontWeight','bold')
text(0.9,0.37,'(std. dev.)','Color','r','FontWeight','bold')
text(0.9,0.24,'(mean)','Color','r','FontWeight','bold')
text(0.9,0.17,'(std. dev.)','Color','r','FontWeight','bold')

% Original Data
text(0.15,0.8,num2str(round(AUC,3)),'Color','b')
text(0.15,0.6,num2str(round(ERRORRATE)),'Color','b')
text(0.15,0.4,num2str(round(PPV,3)),'Color','b')
text(0.15,0.2,num2str(round(PI,3)),'Color','b')

% AM
text(0.35,0.83,num2str(round(mean(AMAUC),3)),'Color','b')
text(0.35,0.63,num2str(round(mean(AME))),'Color','b')
text(0.35,0.43,num2str(round(mean(AMPPV),3)),'Color','b')
text(0.35,0.23,num2str(round(mean(AMPI),3)),'Color','b')
text(0.35,0.77,num2str(round(std(AMAUC),3)),'Color','b')
text(0.35,0.57,num2str(round(std(AME))),'Color','b')
text(0.35,0.37,num2str(round(std(AMPPV),3)),'Color','b')
text(0.35,0.17,num2str(round(std(AMPI),3)),'Color','b')

% MX
text(0.55,0.83,num2str(round(mean(MXAUC),3)),'Color','b')
text(0.55,0.63,num2str(round(mean(MXE))),'Color','b')
text(0.55,0.43,num2str(round(mean(MXPPV),3)),'Color','b')
text(0.55,0.23,num2str(round(mean(MXPI),3)),'Color','b')
text(0.55,0.77,num2str(round(std(MXAUC),3)),'Color','b')
text(0.55,0.57,num2str(round(std(MXE))),'Color','b')
text(0.55,0.37,num2str(round(std(MXPPV),3)),'Color','b')
text(0.55,0.17,num2str(round(std(MXPI),3)),'Color','b')


% GM
text(0.75,0.83,num2str(round(mean(GMAUC),3)),'Color','b')
text(0.75,0.63,num2str(round(mean(GME))),'Color','b')
text(0.75,0.43,num2str(round(mean(GMPPV),3)),'Color','b')
text(0.75,0.23,num2str(round(mean(GMPI),3)),'Color','b')
text(0.75,0.77,num2str(round(std(GMAUC),3)),'Color','b')
text(0.75,0.57,num2str(round(std(GME))),'Color','b')
text(0.75,0.37,num2str(round(std(GMPPV),3)),'Color','b')
text(0.75,0.17,num2str(round(std(GMPI),3)),'Color','b')

annotation('line',[0.25,0.78],[0.65,0.65],'linewidth',2)
annotation('line',[0.25,0.78],[0.5,0.5],'linewidth',2)
annotation('line',[0.25,0.78],[0.34,0.34],'linewidth',2)

