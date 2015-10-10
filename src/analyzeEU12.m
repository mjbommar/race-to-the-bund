clear all;

EU12 = ascii2fts('../data/EU/EU12_crop.txt', 0, 1);

% Determine the list of years in the dataset
yearList = unique(year(EU12.dates));


%% Determine how many yields are available for each year
yearDataAvailable = zeros(length(yearList), 1);
for i = 1:length(yearList)
    cyear = yearList(i);
    cDates = find(year(EU12.dates) == cyear);
    disp(sprintf('%d,%d', cyear, length(cDates)));
    yearDataAvailable(i) = length(cDates);
end

figure;
scatter(yearList, yearDataAvailable);
title('Availability of ten-year yields in the EU-12, 1958-2010');
xlabel('Year');
ylabel('Number of yields available');
print(sprintf('../paper/fig_data_eu12.pdf'), '-dpdf', '-r200');

%% Analyze the spread of the yield data
% Use the German ten-year as the benchmark rate
benchIndex = 3;

matEU12 = fts2mat(EU12);
benchData = matEU12(:, benchIndex);
matEU12 = matEU12(:, setdiff(1:12, [benchIndex]));
spreadData = matEU12 - repmat(benchData, 1, 11);
meanSpreadData = mean(spreadData, 2);
spreadFTS = fints(EU12.dates, spreadData);
meanSpreadFTS = fints(EU12.dates, meanSpreadData);
stdSpreadFTS = fints(EU12.dates, std(spreadData')');

figure;
plot(spreadFTS);
legend({'Austria', 'Belgium', 'Denmark', 'Spain', 'France', 'UK', 'Ireland', 'Italy', 'Netherlands', 'Portugal', 'Sweden'}, 'location', 'best');
title('Spread of EU-12 countries against German ten-year yield');
ylabel('Ten-year yield spread');
xlabel('Date');

figure;
plot(meanSpreadFTS);
title('Mean spread of EU-12 countries against German ten-year yield');
ylabel('Mean ten-year yield spread');
xlabel('Date');

figure;
plot(stdSpreadFTS);
title('Standard deviation of spread of EU-12 countries against German ten-year yield');
ylabel('Std. dev. of ten-year yield spread');
xlabel('Date');

%% Analyze the comovement of the yield data
yearCorrelationEig = zeros(length(yearList), 1);
yearCorrelationEVecStd = zeros(length(yearList), 1);
for i = 1:length(yearList)
    cyear = yearList(i);
    cDates = find(year(EU12.dates) == cyear);
    yieldData = fts2mat(EU12(cDates));
    C = corrcoef(yieldData);
    goodAssets = find(sum(isnan(C)) ~= 12);
    [v,d] = eig(C(goodAssets, goodAssets));
    d = diag(d);
    yearCorrelationEig(i) = max(d);
    yearCorrelationEVecStd(i) = std(v(:,1));
end

% Plot the maximum eigenvalue of the matrix
figure;
scatter(yearList, yearCorrelationEig);
title('Maximum eigenvalue of yield correlation within EU-12, 1958-2010');
xlabel('Year');
ylabel('\lambda_{max}');
print(sprintf('../paper/fig_maxeig_eu12.pdf'), '-dpdf', '-r200');

% Plot the standard deviation of the leading eigenvector
figure;
scatter(yearList, yearCorrelationEVecStd);
title('Standard deviation of the leading eigenvector of yield correlation within EU-12, 1958-2010');
xlabel('Year');
ylabel('\sigma');
print(sprintf('../paper/fig_maxeigstd_eu12.pdf'), '-dpdf', '-r200');

%% Analyze the difference comovement of the yield data
yearCorrelationEig = zeros(length(yearList), 1);
yearCorrelationEVecStd = zeros(length(yearList), 1);

for i = 1:length(yearList)
    cyear = yearList(i);
    cDates = find(year(EU12.dates) == cyear);
    yieldData = diff(fts2mat(EU12(cDates)));
    C = corrcoef(yieldData);
    goodAssets = find(sum(isnan(C)) ~= 12);
    [v,d] = eig(C(goodAssets, goodAssets));
    d = diag(d);
    yearCorrelationEig(i) = max(d);
    yearCorrelationEVecStd(i) = std(v(:,1));
end

% Plot the maximum eigenvalue of the matrix
figure;
scatter(yearList, yearCorrelationEig);
title('Maximum eigenvalue of yield change correlation within EU-12, 1958-2010');
xlabel('Year');
ylabel('\lambda_{max}');
print(sprintf('../paper/fig_diff_maxeig_eu12.pdf'), '-dpdf', '-r200');

% Plot the standard deviation of the leading eigenvector
figure;
scatter(yearList, yearCorrelationEVecStd);
title('Standard deviation of the leading eigenvector of yield change correlation within EU-12, 1958-2010');
xlabel('Year');
ylabel('\sigma');
print(sprintf('../paper/fig_diff_maxeigstd_eu12.pdf'), '-dpdf', '-r200');



%% Discussion section tables
yearRanges = [1958 1966;
    1967 1973;
    1974 1978;
    1979 1985;
    1986 1992;
    1993 1998
    1999 2008;
    2009 2010];

for i = 1:length(yearRanges(:,1))
    startYear = yearRanges(i, 1);
    endYear = yearRanges(i, 2);
    dateIndex = find((year(EU12.dates) >= startYear) .* (year(EU12.dates) <= endYear));
    
    matEU12 = fts2mat(EU12(dateIndex));
    benchData = matEU12(:, benchIndex);
    matEU12 = matEU12(:, setdiff(1:12, [benchIndex]));
    spreadData = abs(matEU12 - repmat(benchData, 1, 11));
    
    disp(sprintf('%d-%d & %0.3f & %0.3f\\\\\\hline', startYear, endYear, mean(mean(spreadData, 2)), mean(std(spreadData')')));
end

for i = 1:length(yearRanges(:,1))
    startYear = yearRanges(i, 1);
    endYear = yearRanges(i, 2);
    dateIndex = find((year(EU12.dates) >= startYear) .* (year(EU12.dates) <= endYear));
    
    
    yieldData = fts2mat(EU12(dateIndex));
    C = corrcoef(yieldData);
    goodAssets = find(sum(isnan(C)) ~= 12);
    [v,d] = eig(C(goodAssets, goodAssets));
    d = diag(d);
    maxEig = max(d);
    
    yieldDiffData = diff(fts2mat(EU12(dateIndex)));
    C = corrcoef(yieldDiffData);
    goodAssets = find(sum(isnan(C)) ~= 12);
    [v,d] = eig(C(goodAssets, goodAssets));
    d = diag(d);
    maxDiffEig = max(d);
    
    disp(sprintf('%d-%d & %0.3f & %0.3f\\\\\\hline', startYear, endYear, maxEig / 12, maxDiffEig / 12));
end

