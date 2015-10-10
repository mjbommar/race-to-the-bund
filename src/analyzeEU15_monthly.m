clear all;

N = 15;
EU15 = ascii2fts('../data/EU/EU15_crop.txt', 0, 1);
matEU15 = fts2mat(EU15) / 100.0;

% Determine the list of years in the dataset
yearList = unique(year(EU15.dates));
benchmarkIndex = 3;

monthDate = zeros(length(yearList) * 12, 1);
dataPoints = zeros(length(yearList) * 12, 1);
spreadMean = zeros(length(yearList) * 12, 1);
spreadStd = zeros(length(yearList) * 12, 1);
spreadBenchmark = matEU15(:, benchmarkIndex);
spreadData = matEU15(:, setdiff(1:N, [benchmarkIndex])) - repmat(spreadBenchmark, 1, N-1);

eigYieldValues = nan(length(yearList) * 12, 1);
eigYieldVectors = nan(length(yearList) * 12, N);
eigDiffYieldValues = nan(length(yearList) * 12, 1);
eigDiffYieldVectors = nan(length(yearList) * 12, N);

for i = 1:length(yearList)
    % Determine the current year
    currentYear = yearList(i);
    
    % Iterate over all months in the year
    for currentMonth = 1:12
        % Store vector index
        vIndex = (i-1) * 12 + currentMonth;
        
        % Calculate the dates that match this year-month constraint
        dateIndex = find((year(EU15.dates) == currentYear) .* (month(EU15.dates) == currentMonth));
        monthStr = sprintf('%02d-%04d', currentMonth, currentYear); 
        monthDate(vIndex) = datenum(monthStr, 'mm-yyyy');
        dataPoints(vIndex) = length(dateIndex);
        
        % Skip months with no data
        if dataPoints(vIndex) == 0
            continue
        end
        
        % Now calculate the yield and yield difference correlation 
        corrYield = corrcoef(log(matEU15(dateIndex, :)));
        goodIndex =  find(sum(isnan(corrYield)) ~= N);
        [V, D] = eig(corrYield(goodIndex, goodIndex));
        eigYieldValues(vIndex) = max(diag(D)) / length(goodIndex);
        eigYieldVectors(vIndex, goodIndex) = V(:,1)';
        
        corrDiffYield = corrcoef(diff(log(matEU15(dateIndex, :))));
        goodIndex =  find(sum(isnan(corrDiffYield)) ~= N);
        [V, D] = eig(corrDiffYield(goodIndex, goodIndex));
        eigDiffYieldValues(vIndex) = max(diag(D)) / length(goodIndex);
        eigDiffYieldVectors(vIndex, goodIndex) = V(:,1)';
    end
end
        