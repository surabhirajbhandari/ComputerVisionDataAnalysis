close all;
clear;
 
[status,sheets] = xlsfinfo('subu.xlsx'); 
[A,names,raw] =xlsread('subu.xlsx',1); 

%% Separate the dataset into Absence (Abs) and Presence (Prs) groups
A; 
Abs=A(1:70);
N0=length(Abs);
Prs=A(71:end); 
N1=length(Prs);

%% Summary of X^2 tests
sorted = sort(A);
xx=0:0.25:max(A);
LN=length(xx);

%% Target Absent
xi=0:size(Abs); 
empirical = ksdensity(sorted,xx);
parWa = fitdist(Abs, 'Weibull');
parNa = fitdist(Abs, 'Nakagami');
parGa = fitdist(Abs, 'Gamma');
parRa = fitdist(Abs, 'Rician');
parLo = fitdist(Abs, 'Lognormal');
parRy = fitdist(Abs, 'Rayleigh');

%% Perform Chi-Square Goodness-of-Fit Tests for "Target Absent" data
[hWa, pWa, statsWa] = chi2gof(Abs, 'CDF',parWa, 'nBins', 7, 'Emin', 1);
X_RW = statsWa.chi2stat / statsWa.df;
[hNa, pNa, statsNa] = chi2gof(Abs, 'CDF',parNa, 'nBins', 7, 'Emin', 1);
X_Gn = statsNa.chi2stat / statsNa.df;
[hGa, pGa, statsGa] = chi2gof(Abs, 'CDF',parGa, 'nBins', 7, 'Emin', 1);
X_Rg = statsGa.chi2stat / statsGa.df;
[hRa, pRa, statsRa] = chi2gof(Abs, 'CDF',parRa, 'nBins', 8, 'Emin', 1);
X_Rr = statsRa.chi2stat / statsRa.df;
[hLo, pLo, statsLo] = chi2gof(Abs, 'CDF', parLo, 'nBins', 7, 'Emin', 1);
X_Lo = statsLo.chi2stat / statsLo.df;
[hRy, pRy, statsRy] = chi2gof(Abs, 'CDF', parRy, 'nBins', 7, 'Emin', 1);
X_Ry = statsRy.chi2stat / statsRy.df

%% Create table summarizing chi-square test results
density = ["Gamma"; "Weibull"; "Nakagami"; "Rician"; "Lognormal"; "Rayleigh"];
h = [hGa; hWa; hNa; hRa; hLo; hRy];
DoF = [statsGa.df; statsWa.df; statsNa.df; statsRa.df; statsLo.df; statsRy.df];
chi_squared_stat = [statsGa.chi2stat; statsWa.chi2stat; statsNa.chi2stat; statsRa.chi2stat; statsLo.chi2stat; statsRy.chi2stat];
p_valuea = [pGa; pWa; pNa; pRa; pLo; pRy];

T_1 = table(h, DoF, chi_squared_stat, p_valuea, 'RowNames', density);
disp("Best fit " + density(find(p_valuea == max(p_valuea))))
disp("Best fit " + density(find(p_valuea == max(p_valuea))))

%% Repeat the same analysis for "Target Present" data
empirical = ksdensity(sorted,xx);
parW = fitdist(Prs, 'Weibull');
parN = fitdist(Prs, 'Nakagami');
parG = fitdist(Prs, 'Gamma');
parR1 = fitdist(Prs, 'Rician');
parL = fitdist(Prs, 'Lognormal');
parRy1 = fitdist(Prs, 'Rayleigh');

[hW, pW, statsW] = chi2gof(Prs, 'CDF',parW, 'nBins', 7, 'Emin', 1);
X_RW = statsW.chi2stat / statsW.df;
[hN, pN, statsN] = chi2gof(Prs, 'CDF',parN, 'nBins', 7, 'Emin', 1);
X_Gn = statsN.chi2stat / statsN.df;
[hG, pG, statsG] = chi2gof(Prs, 'CDF',parG, 'nBins', 7, 'Emin', 1);
X_Rg = statsG.chi2stat / statsG.df;
[hR, pR, statsR] = chi2gof(Prs, 'CDF',parR1, 'nBins', 8, 'Emin', 1);
X_Rr = statsR.chi2stat / statsR.df;
[hL, pL, statsL] = chi2gof(Prs, 'CDF', parL, 'nBins', 7, 'Emin', 1);
X_Lo = statsL.chi2stat / statsL.df;
[hRy1, pRy1, statsRy1] = chi2gof(Prs, 'CDF', parRy1, 'nBins', 7, 'Emin', 1);
X_Ry1 = statsRy1.chi2stat / statsRy1.df;


density = ["Gamma"; "Weibull"; "Nakagami"; "Rician"; "Lognormal"; "Rayleigh"];
h = [hG; hW; hN; hR; hL; hRy1];
DoF = [statsG.df; statsW.df; statsN.df; statsR.df; statsL.df; statsRy1.df];
chi_squared_stat = [statsG.chi2stat; statsW.chi2stat; statsN.chi2stat; statsR.chi2stat; statsL.chi2stat; statsRy1.chi2stat];
p_value = [pG; pW; pN; pR; pL; pRy1];
T_2 = table(h, DoF, chi_squared_stat, p_value, 'RowNames', density)
da=density(find(p_valuea == max(p_valuea)))
d=density(find(p_value == max(p_value)))

fprintf("\t\t\t\t\t\tTarget Absent\n")
disp(T_1)
disp("Best fit " + da)

fprintf("\t\t\t\t\t\tTarget Present\n")
disp(T_2)
disp("Best fit " + d)

%% Theoretical Fits
figure(2)
xx=0:0.1:max(A)*1.7 

% Plot histograms of empirical data
hAbs=histogram(Abs,'normalization', 'pdf')
hold on
hAbs.FaceColor="none"
hAbs.LineWidth=1
hPrs=histogram(Prs,'normalization', 'pdf')
hPrs.FaceColor="#3A3B3C"
hPrs.LineWidth=1

% Fit and plot the best theoretical distributions
parAbs=fitdist(Abs,density(find(p_valuea == max(p_valuea))))
FAbs= pdf(density(find(p_valuea == max(p_valuea))),xx,parAbs.ParameterValues(1),parAbs.ParameterValues(2))
plot(xx,FAbs,'r','LineWidth', 1.5)

parPrs=fitdist(Prs,density(find(p_value == max(p_value))))
FPrs= pdf(density(find(p_value == max(p_value))),xx,parPrs.ParameterValues(1),parPrs.ParameterValues(2))
plot(xx,FPrs,'--b','LineWidth', 1.5)

legend('Data (Target Absent)', 'Data(Target Present)', ['Theoretical Fit (Target Absent): '+density(find(p_valuea == max(p_valuea)))+'('+num2str(parAbs.ParameterValues(1))+','+num2str(parPrs.ParameterValues(2))+')'],['Theoretical Fit (Target Present): '+density(find(p_value == max(p_value)))+'('+num2str(parPrs.ParameterValues(1))+', '+num2str(parPrs.ParameterValues(2))+')'])

%% Receiver Operating Characteristic (ROC) Analysis
figure(3)
thr = 0:0.01:100;
dat=[Abs;Prs];
resp =[zeros(N0,1);ones(N1,1)];

% Compute ROC Curve
[pf,pd,T,AUC,opt]=perfcurve(resp,dat,1);

% Compute Area Under Curve (AUC)
AOC = calculateROCarea(resp,dat) 
A1 = AOC/(2-AOC);
A2 = 2 * AOC^2/(1+AOC);
Std = sqrt((AOC*(1-AOC) + (N1-1)*(A1-AOC^2) + (N0-1)*(A2-AOC^2))/(N0*N1));
areaROC=0.5+polyarea(pf,pd)


[parAbs]=fitdist(Abs,da)
[parPrs]=fitdist(Prs,d)

F_t = cdf(da,xx,parAbs.ParameterValues(1), parAbs.ParameterValues(2)); 
F_t2 = cdf(d,xx,parPrs.ParameterValues(1),parPrs.ParameterValues(2));
PF_t = 1 - F_t;
PD_t = 1 - F_t2;
tArea = polyarea(PF_t,PD_t)+0.5
N3 = length(PF_t);
N4 = length(PD_t);
A3 = tArea/(2-tArea);
A4 = 2 * tArea^2/(1+tArea);

Stdtest = sqrt((tArea*(1-tArea) + (N4-1)*(A3-tArea^2) + (N3-1)*(A4-tArea^2))/(N3*N4))

%ROC EDIT % Gamma Distribution from Part
parg0=fitdist(Abs,'gamma')
parg1=fitdist(Prs,'gamma')
aprs=parg1.a
aabs=parg0.a
bprs=parg1.b
babs=parg0.b
PFG=1-cdf('gamma',thr,aabs,babs)
PDG=1-cdf('gamma',thr,aprs,bprs)
AUCG=0.5 + polyarea(PFG,PDG)


hold on
xlabel("Probability of False Alarm")
ylabel("Probability of Detection")
plot(PF_t,PD_t, 'r--', "linewidth", 1.5)
plot(pf,pd, 'k', 'LineWidth', 1.5) % black line
plot(PFG,PDG,'b-.','linewidth',1.5)

boot_AUC = bootstrp(5000,@calculateROCarea,resp,A);
mean_AUC = mean(boot_AUC);
std_AUC = std(boot_AUC);
nB=5000;% number of boot samples 
rng('default')
NTarget = randi(3,70,1);
YTarget = randi(11,60,1);
Data2=[NTarget;YTarget];
Btarea=bootstrp(nB,@calculateROCarea,resp,Data2);
A1=std(Btarea)
B1=mean(Btarea)


plot(0:1,0:1)
legend(['Theoretical fit: AUC = ',num2str(tArea)],...
    ['Empirical: AUC = ',num2str(AUC)], ...
    ['Bigamma fit: AUC = ',num2str(AUCG)], 'Location', 'NorthEast')

Text = {['AUC(empirical) = ',num2str(AUC),', \sigma = ',num2str(Stdtest)],...
    ['AUC(mean) = ',num2str(mean_AUC),', \sigma = ', num2str(std_AUC),' (bootstrapping)'],...
    ['AUC(mean) = ',num2str(B1),', \sigma = ',num2str(A1),' (parametric bootstrapping)']};
text(0.3,0.2,Text)
box on;
hold off


%% Figure 4
figure(4)
tiledlayout(2,1)
nexttile
histogram(Btarea)
xlim([0.5 1])
a1=prctile(Btarea,2.5) % 2.5%
b1=prctile(Btarea,97.5) % 97.5%
title(['Mean Area=',num2str(B1),' [',num2str(a1),', ',num2str(b1),'], \sigma = ',num2str(A1)])
xlabel('AUC values')
ylabel('frequency')
legend('Parametric bootstrap')
legend('location','northwest')


nexttile
hold on
xlim([0.5 1])
histogram(boot_AUC)
AUC_25=prctile(boot_AUC,2.5) % 2.5%
AUC_975=prctile(boot_AUC,97.5) % 97.5%
title(['Mean Area= ',num2str(mean_AUC),' [',num2str(AUC_25),', ',num2str(AUC_975),'], \sigma = ',num2str(std_AUC)])
xlabel('AUC values')
ylabel('frequency')
legend('Non-Parametric bootstrap')
legend('location','northwest')


%% Bootstrapping for AUC Estimation
nB = 5000;  % Number of bootstrap samples

% Nonparametric Bootstrap
boot_AUC = bootstrp(nB, @calculateROCarea, resp, A);
mean_AUC = mean(boot_AUC);
std_AUC = std(boot_AUC);
ci_AUC = prctile(boot_AUC, [2.5 97.5]); % 95% Confidence Interval

% Parametric Bootstrapping
rng('default')
NTarget = randi(3, 70, 1);
YTarget = randi(11, 60, 1);
Data2 = [NTarget; YTarget];
Btarea = bootstrp(nB, @calculateROCarea, resp, Data2);
mean_Btarea = mean(Btarea);
std_Btarea = std(Btarea);
ci_Btarea = prctile(Btarea, [2.5 97.5]);

%% Figure: Nonparametric & Parametric Bootstrapping
figure;
tiledlayout(2,1);  % Two rows, one column

% ------------- Nonparametric Bootstrapping ROC -------------
nexttile;
hold on;
[pf, pd, ~, AUC, ~] = perfcurve(resp, A, 1);
plot(pf, pd, 'r', 'LineWidth', 2);  % Empirical ROC

% Bootstrap ROC CI (Dashed instead of dense scatter)
for i = 1:500
    boot_idx = randi(length(A), length(A), 1);
    [pf_bs, pd_bs, ~, ~] = perfcurve(resp(boot_idx), A(boot_idx), 1);
    plot(pf_bs, pd_bs, 'b:', 'LineWidth', 0.8);
end

plot(0:1, 0:1, 'k', 'LineWidth', 1.2); % Change to black dashed for clarity


% Title and Labels
title('Nonparametric Bootstrapping', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('1-Specificity', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Sensitivity', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

% CI Annotation
text(0.4, 0.9, sprintf('95%% CI [%0.3f, %0.3f]', ci_AUC(1), ci_AUC(2)), ...
    'FontSize', 12, 'FontWeight', 'bold', 'Color', 'black');

% Legend
legend({'Empirical ROC', '95% CI'}, 'FontSize', 10, 'FontWeight', 'bold', ...
    'Location', 'Southeast', 'Box', 'off');
hold off;

% ------------- Parametric Bootstrapping ROC -------------


nexttile;
hold on;
plot(PF_t, PD_t, 'r', 'LineWidth', 2); % Theoretical Best Fit ROC

% Bootstrap ROC CI (Dashed)
for i = 1:500
    boot_idx = randi(length(A), length(A), 1);
    [pf_bs, pd_bs, ~, ~] = perfcurve(resp(boot_idx), A(boot_idx), 1);
    plot(pf_bs, pd_bs, 'b:', 'LineWidth', 0.8);
end

plot(0:1, 0:1, 'k', 'LineWidth', 1.2); % Change to black dashed for clarity


% Title and Labels
title('Parametric Bootstrapping', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('1-Specificity', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Sensitivity', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

% CI Annotation
text(0.4, 0.9, sprintf('95%% CI [%0.3f, %0.3f]', ci_Btarea(1), ci_Btarea(2)), ...
    'FontSize', 12, 'FontWeight', 'bold', 'Color', 'black');

% Legend
legend({'Best Fit ROC', '95% CI'}, 'FontSize', 10, 'FontWeight', 'bold', ...
    'Location', 'Southeast', 'Box', 'off');
hold off;

%% Helper Function for ROC AUC Calculation
function AOC = calculateROCarea(resp, dat)
    [~, ~, ~, AUC, ~] = perfcurve(resp, dat, 1);
    AOC = AUC;
end




