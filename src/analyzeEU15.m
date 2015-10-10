clear all;

EU15 = ascii2fts('../data/EU/EU15_crop.txt', 0, 1);

% Determine the list of years in the dataset
yearList = unique(year(EU15.dates));


%% Determine how many yields are available for each year
yearDataAvailable = zeros(length(yearList), 1);
for i = 1:length(yearList)
    cyear = yearList(i);
    cDates = find(year(EU15.dates) == cyear);
    disp(sprintf('%d,%d', cyear, length(cDates)));
    yearDataAvailable(i) = length(cDates);
end

% Plot the availability of data for the appendix
figure;
scatter(yearList, yearDataAvailable);
title('Availability of ten-year yields in the EU-15, 1987-2010');
xlabel('Year');
ylabel('Number of yields available');
print(sprintf('../paper/fig_data_eu15.pdf'), '-dpdf', '-r200');


%% Analyze the spread of the yield data
% Use the German ten-year as the benchmark rate
benchIndex = 3;

matEU15 = fts2mat(EU15);
benchData = matEU15(:, benchIndex);
matEU15 = matEU15(:, setdiff(1:15, [benchIndex]));
spreadData = matEU15 - repmat(benchData, 1, 14);
meanSpreadData = mean(spreadData, 2);
spreadFTS = fints(EU15.dates, spreadData);
meanSpreadFTS = fints(EU15.dates, meanSpreadData);
stdSpreadFTS = fints(EU15.dates, std(spreadData')');

figure;
plot(spreadFTS);
legend({'Austria', 'Belgium', 'Denmark', 'Spain', 'Finland', 'France', 'UK', 'Greece', 'Ireland', 'Italy', 'Luxembourg', 'Netherlands', 'Portugal', 'Sweden'}, 'location', 'best');
title('Spread of EU-15 countries against German ten-year yield');
ylabel('Ten-year yield spread');
xlabel('Date');

figure;
plot(meanSpreadFTS);
title('Mean spread of EU-15 countries against German ten-year yield');
ylabel('Mean ten-year yield spread');
xlabel('Date');

figure;
plot(stdSpreadFTS);
title('Standard deviation of spread of EU-15 countries against German ten-year yield');
ylabel('Std. dev. of ten-year yield spread');
xlabel('Date');

%% Analyze the comovement of the yield data
yearCorrelationEig = zeros(length(yearList), 1);
yearCorrelationEVecStd = zeros(length(yearList), 1);
for i = 1:length(yearList)
    cyear = yearList(i);
    cDates = find(year(EU15.dates) == cyear);
    yieldData = fts2mat(EU15(cDates));
    C = corrcoef(yieldData);
    goodAssets = find(sum(isnan(C)) ~= 15);
    [v,d] = eig(C(goodAssets, goodAssets));
    d = diag(d);
    yearCorrelationEig(i) = max(d);
    yearCorrelationEVecStd(i) = std(v(:,1));
end

% Plot the maximum eigenvalue of the matrix
figure;
scatter(yearList, yearCorrelationEig);
title('Maximum eigenvalue of yield correlation within EU-15, 1987-2010');
xlabel('Year');
ylabel('\lambda_{max}');
print(sprintf('../paper/fig_maxeig_eu15.pdf'), '-dpdf', '-r200');

% Plot the standard deviation of the leading eigenvector
figure;
scatter(yearList, yearCorrelationEVecStd);
title('Standard deviation of the leading eigenvector of yield correlation within EU-15, 1987-2010');
xlabel('Year');
ylabel('\sigma');
print(sprintf('../paper/fig_maxeigstd_eu15.pdf'), '-dpdf', '-r200');

%% Analyze the difference comovement of the yield data
yearCorrelationEig = zeros(length(yearList), 1);
yearCorrelationEVecStd = zeros(length(yearList), 1);

for i = 1:length(yearList)
    cyear = yearList(i);
    cDates = find(year(EU15.dates) == cyear);
    yieldData = diff(fts2mat(EU15(cDates)));
    C = corrcoef(yieldData);
    goodAssets = find(sum(isnan(C)) ~= 15);
    [v,d] = eig(C(goodAssets, goodAssets));
    d = diag(d);
    yearCorrelationEig(i) = max(d);
    yearCorrelationEVecStd(i) = std(v(:,1));
end

% Plot the maximum eigenvalue of the matrix
figure;
scatter(yearList, yearCorrelationEig);
title('Maximum eigenvalue of yield change correlation within EU-15, 1987-2010');
xlabel('Year');
ylabel('\lambda_{max}');
print(sprintf('../paper/fig_diff_maxeig_eu15.pdf'), '-dpdf', '-r200');

% Plot the standard deviation of the leading eigenvector
figure;
scatter(yearList, yearCorrelationEVecStd);
title('Standard deviation of the leading eigenvector of yield change correlation within EU-15, 1987-2010');
xlabel('Year');
ylabel('\sigma');
print(sprintf('../paper/fig_diff_maxeigstd_eu15.pdf'), '-dpdf', '-r200');