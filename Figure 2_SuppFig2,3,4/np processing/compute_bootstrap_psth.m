function NPdata = compute_bootstrap_psth(NPdata, nBoot, conditionField)
% Compute bootstrapped mean PSTH and SEM from raster data in NPdata

for animalIdx = 1:numel(NPdata)
    units = NPdata(animalIdx).cells.units;
    for unitIdx = 1:numel(units)
        % Get raster data directly
        rasterData = units(unitIdx).(conditionField).Raster;  % [time bins x trials]
        
        [nTimeBins, nTrials] = size(rasterData);
        bootPSTH = zeros(nTimeBins, nBoot);
        
        for b = 1:nBoot
            % Sample trials with replacement
            resampledIdx = randi(nTrials, [1, nTrials]);
            resampledData = rasterData(:, resampledIdx);
            
            % Average across sampled trials (mean at each time bin)
            bootPSTH(:, b) = mean(resampledData, 2, 'omitnan');
        end
        
        % Compute mean and SEM across bootstrapped PSTHs
        meanPSTH = mean(bootPSTH, 2);
        semPSTH = std(bootPSTH, 0, 2);
        
        % Save back into structure
        NPdata(animalIdx).cells.units(unitIdx).(conditionField).meanPSTH = meanPSTH;
        NPdata(animalIdx).cells.units(unitIdx).(conditionField).semPSTH = semPSTH;
    end
end
end
