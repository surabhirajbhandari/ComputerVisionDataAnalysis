% Load data from Excel
data = xlsread('Subu.xlsx');
H0 = xlsread('Subu.xlsx', 1, 'A1:A70'); % Target Absent
H1 = xlsread('Subu.xlsx', 1, 'A71:A130'); % Target Present

% Set axis limits
xlim([0,1.2*max(data)]);
hold on

% Plot histograms and kernel density estimates (KDEs)
histogram(H0,'Normalization','pdf');
[f,xi] = ksdensity(H0);
plot(xi,f, 'b', 'LineWidth', 1.5);
[mf,mxi] = max(f);
plot(xi(find(f == max(f))),mf,'bo', 'MarkerSize', 8, 'LineWidth', 2);

histogram(H1,'Normalization','pdf');
[f2,xj] = ksdensity(H1);
plot(xj,f2, 'r--', 'LineWidth', 1.5);
[mf2,mxj] = max(f2);
plot(xj(find(f2 == max(f2))),mf2,'bs', 'MarkerSize', 8, 'LineWidth', 2);


[inter_x, inter_y] = polyxpoly(xi, f, xj, f2);

if isempty(inter_x)
    error('No intersection found between the PDFs. Check the data distribution.');
end

% Select the second intersection as the threshold dynamically
ThreshInt = inter_x(end); 

legend('H0', 'fit(H0)', 'max(H0)', 'H1', 'fit(H1)', 'max(H1)')

xlabel('Values');
ylabel('Estimated Density');

%----------------- HW3: CDF PLOT ------------------
hold off
figure(2)
ecdf(H0);
hold on
ecdf(H1);

% Generate x-axis values
xax = linspace(0, max(data), 1000);

% Compute CDFs using ksdensity
FH0 = ksdensity(H0, xax, 'Function', 'CDF');
FH1 = ksdensity(H1, xax, 'Function', 'CDF');

% Plot smooth CDFs
plot(xax, FH0, 'b', 'LineWidth', 1.5);
plot(xax, FH1, 'r', 'LineWidth', 1.5);

% **FIX:** Use interpolation instead of `find()`
PF_CDF = 1 - interp1(xax, FH0, ThreshInt, 'linear', 'extrap');
PM_CDF = interp1(xax, FH1, ThreshInt, 'linear', 'extrap');

% Mark points for P_F and P_M
plot(ThreshInt, PF_CDF, 'rs', 'MarkerSize', 10, 'LineWidth', 2);
plot(ThreshInt, PM_CDF, 'b*', 'MarkerSize', 10, 'LineWidth', 2);

% Labeling and legends
legend('ecdf(c|H0)','ecdf(c|H0)','F_v(v|H0)','F_v(v|H1)','1-P_F','P_M')

xlabel('Values (v)');
ylabel('Estimated CDF');




% Display probabilities with proper formatting
text(ThreshInt + 0.25, 0.15, ['P_F = ', num2str(1 - PF_CDF, '%.4f'), ' (count), ', num2str(PF_CDF, '%.4f'), ' (CDF)'], 'FontWeight', 'bold');
text(ThreshInt + 0.25, 0.05, ['P_M = ', num2str(PM_CDF, '%.4f'), ' (count), ', num2str(PM_CDF, '%.4f'), ' (CDF)'], 'FontWeight', 'bold');

grid on;
hold off;

