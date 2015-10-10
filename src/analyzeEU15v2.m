clear all;

N = 12;
EU15 = ascii2fts('../data/EU/EU12_crop.txt', 0, 1);
matEU15 = fts2mat(EU15) / 100.0;

% Determine the list of years in the dataset
yearList = unique(year(EU15.dates));
benchmarkIndex = 3;

monthDate = zeros(length(yearList), 1);
dataPoints = zeros(length(yearList), 1);
spreadBenchmark = matEU15(:, benchmarkIndex);
spreadData = matEU15(:, setdiff(1:N, [benchmarkIndex])) - repmat(spreadBenchmark, 1, N-1);
spreadMean = mean(spreadData, 2);
spreadStd = std(spreadData')';

eigYieldValues = nan(length(yearList), 1);
eigYieldVectors = nan(length(yearList), N);
eigDiffYieldValues = nan(length(yearList), 1);
eigDiffYieldVectors = nan(length(yearList), N);

for i = 1:length(yearList)
    % Determine the current year
    currentYear = yearList(i);
    
    % Calculate the dates that match this year-month constraint
    dateIndex = find((year(EU15.dates) == currentYear));
    dataPoints(i) = length(dateIndex);
        
    % Skip months with no data
    if dataPoints(i) == 0
        continue
    end
        
    % Now calculate the yield and yield difference correlation 
    corrYield = corrcoef(log(matEU15(dateIndex, :)));
    goodIndex =  find(sum(isnan(corrYield)) ~= N);
    [V, D] = eig(corrYield(goodIndex, goodIndex));
    eigYieldValues(i) = max(diag(D)) / length(goodIndex);
    eigYieldVectors(i, goodIndex) = V(:,1)';

    corrDiffYield = corrcoef(diff(log(matEU15(dateIndex, :))));
    goodIndex =  find(sum(isnan(corrDiffYield)) ~= N);
    [V, D] = eig(corrDiffYield(goodIndex, goodIndex));
    eigDiffYieldValues(i) = max(diag(D)) / length(goodIndex);
    eigDiffYieldVectors(i, goodIndex) = V(:,1)';
end
        