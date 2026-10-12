

%% ========================================================================
% DECELERATION-AWARE ENDPOINT PREDICTION FOR CSM REACHES
%
% For every T0-T3 reach:
%
%   1. Estimate velocity and acceleration immediately before CSM onset.
%   2. Predict when outward velocity would reach zero without a CSM.
%   3. Predict the x-y stopping position at that time.
%   4. Define the actual endpoint as maximum outward x position.
%   5. Compare predicted and actual endpoint error relative to the target.
%
% DATA ORGANIZATION
%
%   dataStruct.(animalName).ST0_CSM.stimvel
%   dataStruct.(animalName).ST1_CSM.stimvel
%   dataStruct.(animalName).ST2_CSM.stimvel
%   dataStruct.(animalName).ST3_CSM.stimvel
%
% stimvel dimensions:
%   samples x variables x reaches
%
% Columns:
%   1   = time
%   2:4 = position [x y z], in decimeters
%   5:7 = velocity [vx vy vz]
%
% CSM onset is aligned at row 21.
% ========================================================================

clearvars -except ExcitRR_Sham ExcitRR_Ctrl
close all
clc


%% ========================================================================
% SETTINGS
% ========================================================================

dataStruct = ExcitRR_Sham;

subtypes = {'T0','T1','T2','T3'};

timeColumn = 1;
positionColumns = 2:4;
velocityColumns = 5:7;

% Analyze x-y plane.
analysisDimensions = [1 2];

% Outward direction is x.
outwardDimension = 1;

% CSM onset.
CSM_idx = 21;

% Window used to estimate pre-CSM velocity and acceleration.
%
% Eight samples gives a more stable acceleration estimate than five.
preFitWindow = (CSM_idx-8):(CSM_idx-1);

% Number of final pre-CSM samples used to estimate velocity at onset.
velocityAtOnsetWindow = (CSM_idx-3):(CSM_idx-1);

% Time units in column 1:
%   'seconds'
%   'milliseconds'
%   'auto'
timeUnits = 'auto';

% Only accept predicted stopping times within this range.
%
% These values are in seconds.
minimumPredictedStopTime = 0.005;
maximumPredictedStopTime = 0.250;

% Require outward velocity before CSM to be positive.
minimumOutwardVelocity = 0;

% Require outward acceleration to be negative.
%
% This is necessary for a future maximum-outward stopping point.
requireOutwardDeceleration = true;

minimumReachesPerAnimal = 5;

% Plot in centimeters.
positionScale = 10;
positionUnit = 'cm';

% Enter the measured target location [x y z] in decimeters.
%
% Leaving this as NaN estimates the target from the median actual endpoint
% across every subtype and animal.
targetPosition = [NaN NaN NaN];

% Ellipse settings.
ellipseScale = 2;
numberEllipsePoints = 300;

rng(1);


%% ========================================================================
% IDENTIFY ANIMALS
% ========================================================================

animalNames = fieldnames(dataStruct);
numberAnimals = numel(animalNames);
numberSubtypes = numel(subtypes);


%% ========================================================================
% ESTIMATE TARGET IF NECESSARY
% ========================================================================

if any(~isfinite(targetPosition))

    candidateEndpoints = [];

    for subtypeIdx = 1:numberSubtypes

        subtype = subtypes{subtypeIdx};
        csmField = ['S' subtype '_CSM'];

        for animalIdx = 1:numberAnimals

            animalName = animalNames{animalIdx};

            if ~isfield(dataStruct.(animalName),csmField)
                continue
            end

            subtypeStruct = ...
                dataStruct.(animalName).(csmField);

            if ~isfield(subtypeStruct,'stimvel') || ...
                    isempty(subtypeStruct.stimvel)
                continue
            end

            stimvel = subtypeStruct.stimvel;

            if ~isnumeric(stimvel) || ndims(stimvel) ~= 3
                continue
            end

            for reachIdx = 1:size(stimvel,3)

                reach = stimvel(:,:,reachIdx);

                if size(reach,2) < max(positionColumns)
                    continue
                end

                position3D = reach(:,positionColumns);
                outwardPosition = ...
                    position3D(:,outwardDimension);

                finiteOutward = isfinite(outwardPosition);

                if ~any(finiteOutward)
                    continue
                end

                outwardForMaximum = outwardPosition;
                outwardForMaximum(~finiteOutward) = -Inf;

                [~,endpointIdx] = max(outwardForMaximum);

                endpoint = position3D(endpointIdx,:);

                if all(isfinite(endpoint))
                    candidateEndpoints(end+1,:) = endpoint;
                end
            end
        end
    end

    if isempty(candidateEndpoints)
        error('No endpoints were available for target estimation.');
    end

    targetPosition = ...
        median(candidateEndpoints,1,'omitnan');

    fprintf('\nTarget estimated from pooled actual endpoints.\n');

else

    fprintf('\nUsing manually specified target position.\n');

end

fprintf('Target position [x y z], dm:\n');
disp(targetPosition);


%% ========================================================================
% RUN ANALYSIS
% ========================================================================

DecelerationEndpointResults = struct;

for subtypeIdx = 1:numberSubtypes

    subtype = subtypes{subtypeIdx};
    csmField = ['S' subtype '_CSM'];

    fprintf('\n============================================================\n');
    fprintf('ANALYZING %s\n',subtype);
    fprintf('============================================================\n');

    %% Animal-level arrays

    animalRawReachCount = zeros(numberAnimals,1);
    animalValidReachCount = zeros(numberAnimals,1);

    animalMeanPredictedError = nan(numberAnimals,1);
    animalMeanActualError = nan(numberAnimals,1);
    animalMeanImprovement = nan(numberAnimals,1);
    animalMeanTargetProjection = nan(numberAnimals,1);
    animalFractionImproved = nan(numberAnimals,1);

    %% Exclusion audit

    excludedNoField = zeros(numberAnimals,1);
    excludedShortReach = zeros(numberAnimals,1);
    excludedInvalidData = zeros(numberAnimals,1);
    excludedNotMovingOutward = zeros(numberAnimals,1);
    excludedNotDecelerating = zeros(numberAnimals,1);
    excludedStopTime = zeros(numberAnimals,1);
    excludedEndpointBeforeCSM = zeros(numberAnimals,1);

    %% Reach-level pooled arrays

    allPredictedEndpoint = [];
    allActualEndpoint = [];
    allOnsetPosition = [];

    allPredictedError = [];
    allActualError = [];
    allImprovement = [];
    allTargetProjection = [];

    allPredictedStopTime = [];
    allActualEndpointTime = [];

    allVelocityAtOnset = [];
    allAcceleration = [];

    allAnimalID = [];
    allOriginalReachIndex = [];

    target2D = targetPosition(analysisDimensions);


    %% --------------------------------------------------------------------
    % Animal loop
    % ---------------------------------------------------------------------

    for animalIdx = 1:numberAnimals

        animalName = animalNames{animalIdx};

        if ~isfield(dataStruct.(animalName),csmField)

            excludedNoField(animalIdx) = 1;
            continue

        end

        subtypeStruct = dataStruct.(animalName).(csmField);

        if ~isfield(subtypeStruct,'stimvel') || ...
                isempty(subtypeStruct.stimvel)
            continue
        end

        stimvel = subtypeStruct.stimvel;

        if ~isnumeric(stimvel) || ndims(stimvel) ~= 3
            continue
        end

        numberReaches = size(stimvel,3);
        animalRawReachCount(animalIdx) = numberReaches;

        predictedEndpoint = nan(numberReaches,2);
        actualEndpoint = nan(numberReaches,2);
        onsetPosition = nan(numberReaches,2);

        predictedError = nan(numberReaches,1);
        actualError = nan(numberReaches,1);
        improvement = nan(numberReaches,1);
        targetProjection = nan(numberReaches,1);

        predictedStopTime = nan(numberReaches,1);
        actualEndpointTime = nan(numberReaches,1);

        velocityAtOnset = nan(numberReaches,2);
        accelerationAtOnset = nan(numberReaches,2);


        %% ---------------------------------------------------------------
        % Reach loop
        % ---------------------------------------------------------------

        for reachIdx = 1:numberReaches

            reach = stimvel(:,:,reachIdx);

            if isempty(reach) || ...
                    size(reach,1) < CSM_idx || ...
                    size(reach,2) < max(velocityColumns) || ...
                    min(preFitWindow) < 1 || ...
                    max(preFitWindow) > size(reach,1)

                excludedShortReach(animalIdx) = ...
                    excludedShortReach(animalIdx) + 1;

                continue
            end

            timeVector = reach(:,timeColumn);
            position3D = reach(:,positionColumns);
            velocity3D = reach(:,velocityColumns);

            position2D = position3D(:,analysisDimensions);
            velocity2D = velocity3D(:,analysisDimensions);

            %% Convert full time vector to seconds

            timeSeconds = convertTimeVectorToSecondsLocal( ...
                timeVector,timeUnits);

            if any(~isfinite(timeSeconds(preFitWindow))) || ...
                    any(any(~isfinite( ...
                    velocity2D(preFitWindow,:)))) || ...
                    any(~isfinite(position2D(CSM_idx,:)))

                excludedInvalidData(animalIdx) = ...
                    excludedInvalidData(animalIdx) + 1;

                continue
            end

            onsetTime = timeSeconds(CSM_idx);
            currentOnsetPosition = position2D(CSM_idx,:);

            %% Estimate velocity at CSM onset

            currentVelocityAtOnset = mean( ...
                velocity2D(velocityAtOnsetWindow,:), ...
                1,'omitnan');

            %% Estimate acceleration by fitting velocity versus time

            fitTimes = ...
                timeSeconds(preFitWindow) - onsetTime;

            currentAcceleration = nan(1,2);

            for dimensionIdx = 1:2

                velocityFit = polyfit( ...
                    fitTimes, ...
                    velocity2D(preFitWindow,dimensionIdx), ...
                    1);

                % Slope of velocity-versus-time is acceleration.
                currentAcceleration(dimensionIdx) = ...
                    velocityFit(1);
            end

            if any(~isfinite(currentVelocityAtOnset)) || ...
                    any(~isfinite(currentAcceleration))

                excludedInvalidData(animalIdx) = ...
                    excludedInvalidData(animalIdx) + 1;

                continue
            end

            outwardVelocity = ...
                currentVelocityAtOnset(outwardDimension);

            outwardAcceleration = ...
                currentAcceleration(outwardDimension);

            %% Must still be moving outward at CSM onset

            if outwardVelocity <= minimumOutwardVelocity

                excludedNotMovingOutward(animalIdx) = ...
                    excludedNotMovingOutward(animalIdx) + 1;

                continue
            end

            %% Must be decelerating outward

            if requireOutwardDeceleration && ...
                    outwardAcceleration >= 0

                excludedNotDecelerating(animalIdx) = ...
                    excludedNotDecelerating(animalIdx) + 1;

                continue
            end

            %% Predicted time until outward velocity reaches zero

            currentPredictedStopTime = ...
                -outwardVelocity / outwardAcceleration;

            if ~isfinite(currentPredictedStopTime) || ...
                    currentPredictedStopTime < ...
                    minimumPredictedStopTime || ...
                    currentPredictedStopTime > ...
                    maximumPredictedStopTime

                excludedStopTime(animalIdx) = ...
                    excludedStopTime(animalIdx) + 1;

                continue
            end

            %% Deceleration-aware predicted endpoint

            currentPredictedEndpoint = ...
                currentOnsetPosition + ...
                currentVelocityAtOnset .* ...
                currentPredictedStopTime + ...
                0.5 .* currentAcceleration .* ...
                currentPredictedStopTime.^2;

            %% Actual endpoint at maximum outward x

            outwardPosition = ...
                position3D(:,outwardDimension);

            finiteOutward = isfinite(outwardPosition);

            if ~any(finiteOutward)

                excludedInvalidData(animalIdx) = ...
                    excludedInvalidData(animalIdx) + 1;

                continue
            end

            outwardForMaximum = outwardPosition;
            outwardForMaximum(~finiteOutward) = -Inf;

            [~,actualEndpointIdx] = max(outwardForMaximum);

            if actualEndpointIdx <= CSM_idx

                excludedEndpointBeforeCSM(animalIdx) = ...
                    excludedEndpointBeforeCSM(animalIdx) + 1;

                continue
            end

            currentActualEndpoint = ...
                position2D(actualEndpointIdx,:);

            currentActualEndpointTime = ...
                timeSeconds(actualEndpointIdx) - onsetTime;

            if any(~isfinite(currentPredictedEndpoint)) || ...
                    any(~isfinite(currentActualEndpoint)) || ...
                    ~isfinite(currentActualEndpointTime)

                excludedInvalidData(animalIdx) = ...
                    excludedInvalidData(animalIdx) + 1;

                continue
            end

            %% Target errors

            currentPredictedError = ...
                norm(currentPredictedEndpoint-target2D);

            currentActualError = ...
                norm(currentActualEndpoint-target2D);

            % Positive means the actual endpoint is closer to the target.
            currentImprovement = ...
                currentPredictedError-currentActualError;

            %% Target-directed component of endpoint shift

            endpointShift = ...
                currentActualEndpoint-currentPredictedEndpoint;

            targetVector = ...
                target2D-currentPredictedEndpoint;

            targetDistance = norm(targetVector);

            if targetDistance > 0

                targetUnitVector = ...
                    targetVector/targetDistance;

                currentTargetProjection = ...
                    dot(endpointShift,targetUnitVector);

            else

                currentTargetProjection = 0;

            end

            %% Store reach

            predictedEndpoint(reachIdx,:) = ...
                currentPredictedEndpoint;

            actualEndpoint(reachIdx,:) = ...
                currentActualEndpoint;

            onsetPosition(reachIdx,:) = ...
                currentOnsetPosition;

            predictedError(reachIdx) = ...
                currentPredictedError;

            actualError(reachIdx) = ...
                currentActualError;

            improvement(reachIdx) = ...
                currentImprovement;

            targetProjection(reachIdx) = ...
                currentTargetProjection;

            predictedStopTime(reachIdx) = ...
                currentPredictedStopTime;

            actualEndpointTime(reachIdx) = ...
                currentActualEndpointTime;

            velocityAtOnset(reachIdx,:) = ...
                currentVelocityAtOnset;

            accelerationAtOnset(reachIdx,:) = ...
                currentAcceleration;
        end


        %% Valid reach mask

        validReach = ...
            all(isfinite(predictedEndpoint),2) & ...
            all(isfinite(actualEndpoint),2) & ...
            isfinite(predictedError) & ...
            isfinite(actualError) & ...
            isfinite(improvement) & ...
            isfinite(targetProjection);

        numberValid = sum(validReach);

        animalValidReachCount(animalIdx) = numberValid;

        fprintf([ ...
            '%s: raw = %d, valid = %d, ' ...
            'not outward = %d, not decelerating = %d, ' ...
            'invalid stop time = %d\n'], ...
            animalName, ...
            numberReaches, ...
            numberValid, ...
            excludedNotMovingOutward(animalIdx), ...
            excludedNotDecelerating(animalIdx), ...
            excludedStopTime(animalIdx));


        %% Animal-level statistics

        if numberValid >= minimumReachesPerAnimal

            animalMeanPredictedError(animalIdx) = ...
                mean(predictedError(validReach),'omitnan');

            animalMeanActualError(animalIdx) = ...
                mean(actualError(validReach),'omitnan');

            animalMeanImprovement(animalIdx) = ...
                mean(improvement(validReach),'omitnan');

            animalMeanTargetProjection(animalIdx) = ...
                mean(targetProjection(validReach),'omitnan');

            animalFractionImproved(animalIdx) = ...
                mean(improvement(validReach)>0);
        end


        %% Pool valid reaches

        originalReachIndices = find(validReach);

        allPredictedEndpoint = [ ...
            allPredictedEndpoint; ...
            predictedEndpoint(validReach,:)];

        allActualEndpoint = [ ...
            allActualEndpoint; ...
            actualEndpoint(validReach,:)];

        allOnsetPosition = [ ...
            allOnsetPosition; ...
            onsetPosition(validReach,:)];

        allPredictedError = [ ...
            allPredictedError; ...
            predictedError(validReach)];

        allActualError = [ ...
            allActualError; ...
            actualError(validReach)];

        allImprovement = [ ...
            allImprovement; ...
            improvement(validReach)];

        allTargetProjection = [ ...
            allTargetProjection; ...
            targetProjection(validReach)];

        allPredictedStopTime = [ ...
            allPredictedStopTime; ...
            predictedStopTime(validReach)];

        allActualEndpointTime = [ ...
            allActualEndpointTime; ...
            actualEndpointTime(validReach)];

        allVelocityAtOnset = [ ...
            allVelocityAtOnset; ...
            velocityAtOnset(validReach,:)];

        allAcceleration = [ ...
            allAcceleration; ...
            accelerationAtOnset(validReach,:)];

        allAnimalID = [ ...
            allAnimalID; ...
            repmat(animalIdx,numberValid,1)];

        allOriginalReachIndex = [ ...
            allOriginalReachIndex; ...
            originalReachIndices];
    end


    %% --------------------------------------------------------------------
    % Animal-level inference
    % ---------------------------------------------------------------------

    validAnimalsError = ...
        isfinite(animalMeanPredictedError) & ...
        isfinite(animalMeanActualError);

    validAnimalsImprovement = ...
        isfinite(animalMeanImprovement);

    validAnimalsProjection = ...
        isfinite(animalMeanTargetProjection);

    validAnimalsFraction = ...
        isfinite(animalFractionImproved);


    if sum(validAnimalsError)>=2

        pPredictedVsActual = signrank( ...
            animalMeanActualError(validAnimalsError), ...
            animalMeanPredictedError(validAnimalsError), ...
            'tail','left');

    else

        pPredictedVsActual = NaN;

    end


    if sum(validAnimalsImprovement)>=2

        pImprovement = signrank( ...
            animalMeanImprovement(validAnimalsImprovement), ...
            0, ...
            'tail','right');

    else

        pImprovement = NaN;

    end


    if sum(validAnimalsProjection)>=2

        pTargetProjection = signrank( ...
            animalMeanTargetProjection(validAnimalsProjection), ...
            0, ...
            'tail','right');

    else

        pTargetProjection = NaN;

    end


    if sum(validAnimalsFraction)>=2

        pFractionImproved = signrank( ...
            animalFractionImproved(validAnimalsFraction), ...
            0.5, ...
            'tail','right');

    else

        pFractionImproved = NaN;

    end


    %% --------------------------------------------------------------------
    % Audit table
    % ---------------------------------------------------------------------

    Animal = animalNames;
    RawReaches = animalRawReachCount;
    ValidReaches = animalValidReachCount;

    NotMovingOutward = excludedNotMovingOutward;
    NotDecelerating = excludedNotDecelerating;
    InvalidStopTime = excludedStopTime;
    EndpointBeforeCSM = excludedEndpointBeforeCSM;
    InvalidData = excludedInvalidData;
    ShortReach = excludedShortReach;

    AuditTable = table( ...
        Animal, ...
        RawReaches, ...
        ValidReaches, ...
        NotMovingOutward, ...
        NotDecelerating, ...
        InvalidStopTime, ...
        EndpointBeforeCSM, ...
        InvalidData, ...
        ShortReach);


    %% --------------------------------------------------------------------
    % Store subtype results
    % ---------------------------------------------------------------------

    R = struct;

    R.subtype = subtype;
    R.rawReachCount = sum(animalRawReachCount);
    R.validReachCount = size(allPredictedEndpoint,1);

    R.animalNames = animalNames;
    R.animalRawReachCount = animalRawReachCount;
    R.animalReachCount = animalValidReachCount;

    R.validAnimalsError = validAnimalsError;
    R.validAnimalsImprovement = validAnimalsImprovement;
    R.validAnimalsProjection = validAnimalsProjection;

    R.animalMeanPredictedError = ...
        animalMeanPredictedError;

    R.animalMeanActualError = ...
        animalMeanActualError;

    R.animalMeanImprovement = ...
        animalMeanImprovement;

    R.animalMeanTargetProjection = ...
        animalMeanTargetProjection;

    R.animalFractionImproved = ...
        animalFractionImproved;

    R.reach.predictedEndpoint = ...
        allPredictedEndpoint;

    R.reach.actualEndpoint = ...
        allActualEndpoint;

    R.reach.onsetPosition = ...
        allOnsetPosition;

    R.reach.predictedError = ...
        allPredictedError;

    R.reach.actualError = ...
        allActualError;

    R.reach.improvement = ...
        allImprovement;

    R.reach.targetProjection = ...
        allTargetProjection;

    R.reach.predictedStopTime = ...
        allPredictedStopTime;

    R.reach.actualEndpointTime = ...
        allActualEndpointTime;

    R.reach.velocityAtOnset = ...
        allVelocityAtOnset;

    R.reach.acceleration = ...
        allAcceleration;

    R.reach.animalID = ...
        allAnimalID;

    R.reach.originalReachIndex = ...
        allOriginalReachIndex;

    R.statistics.pPredictedVsActual = ...
        pPredictedVsActual;

    R.statistics.pImprovement = ...
        pImprovement;

    R.statistics.pTargetProjection = ...
        pTargetProjection;

    R.statistics.pFractionImproved = ...
        pFractionImproved;

    R.auditTable = AuditTable;

    DecelerationEndpointResults.(subtype) = R;
end


%% ========================================================================
% PRINT RESULTS
% ========================================================================

fprintf('\n============================================================\n');
fprintf('DECELERATION-AWARE RESULTS\n');
fprintf('============================================================\n');

for subtypeIdx = 1:numberSubtypes

    subtype = subtypes{subtypeIdx};
    R = DecelerationEndpointResults.(subtype);

    fprintf('\n%s\n',subtype);
    fprintf('Raw reaches: %d\n',R.rawReachCount);
    fprintf('Valid model reaches: %d\n',R.validReachCount);

    fprintf('Mean predicted endpoint error: %.4f dm\n', ...
        mean( ...
        R.animalMeanPredictedError( ...
        R.validAnimalsError), ...
        'omitnan'));

    fprintf('Mean actual endpoint error: %.4f dm\n', ...
        mean( ...
        R.animalMeanActualError( ...
        R.validAnimalsError), ...
        'omitnan'));

    fprintf('Predicted versus actual p = %.6g\n', ...
        R.statistics.pPredictedVsActual);

    fprintf('Mean improvement: %.4f dm, p = %.6g\n', ...
        mean( ...
        R.animalMeanImprovement( ...
        R.validAnimalsImprovement), ...
        'omitnan'), ...
        R.statistics.pImprovement);

    fprintf('Mean target projection: %.4f dm, p = %.6g\n', ...
        mean( ...
        R.animalMeanTargetProjection( ...
        R.validAnimalsProjection), ...
        'omitnan'), ...
        R.statistics.pTargetProjection);
end


%% ========================================================================
% FIGURE 1: PREDICTED VERSUS ACTUAL ENDPOINT ERROR
% ========================================================================

figure( ...
    'Name','Deceleration-aware endpoint error', ...
    'Color','w', ...
    'Position',[100 100 1200 900]);

layoutError = tiledlayout( ...
    2,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

for subtypeIdx = 1:numberSubtypes

    subtype = subtypes{subtypeIdx};
    R = DecelerationEndpointResults.(subtype);

    ax = nexttile;
    hold(ax,'on');

    validAnimals = R.validAnimalsError;

    if sum(validAnimals)>0

        plotPairedValuesLocal( ...
            R.animalMeanPredictedError(validAnimals) .* ...
            positionScale, ...
            R.animalMeanActualError(validAnimals) .* ...
            positionScale);
    end

    xticks(ax,[1 2]);
    xticklabels(ax,{'Predicted','Actual'});

    ylabel(ax,sprintf('Endpoint error (%s)',positionUnit));

    title(ax,sprintf( ...
        '%s: p = %.3g', ...
        subtype, ...
        R.statistics.pPredictedVsActual));

    box(ax,'off');
end




%% ========================================================================
% FIGURE 3: PREDICTED VERSUS ACTUAL ENDPOINT TIME
%
% This validates whether the deceleration model predicts a plausible
% stopping time.
% ========================================================================

figure( ...
    'Name','Predicted versus actual endpoint time', ...
    'Color','w', ...
    'Position',[150 100 1200 900]);

layoutTime = tiledlayout( ...
    2,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

for subtypeIdx = 1:numberSubtypes

    subtype = subtypes{subtypeIdx};
    R = DecelerationEndpointResults.(subtype);

    ax = nexttile;
    hold(ax,'on');

    predictedTimeMS = ...
        R.reach.predictedStopTime .* 1000;

    actualTimeMS = ...
        R.reach.actualEndpointTime .* 1000;

    scatter( ...
        ax, ...
        predictedTimeMS, ...
        actualTimeMS, ...
        20, ...
        'k', ...
        'filled', ...
        'MarkerFaceAlpha',0.25);

    combinedTime = [predictedTimeMS;actualTimeMS];
    combinedTime = combinedTime(isfinite(combinedTime));

    if ~isempty(combinedTime)

        timeLimits = [min(combinedTime),max(combinedTime)];

        plot(ax,timeLimits,timeLimits,'k--','LineWidth',1.5);

        xlim(ax,timeLimits);
        ylim(ax,timeLimits);
    end

    axis(ax,'square');

    xlabel(ax,'Predicted stopping time (ms)');
    ylabel(ax,'Actual endpoint time (ms)');

    title(ax,sprintf('%s, n = %d', ...
        subtype,numel(predictedTimeMS)));

    box(ax,'off');
end


%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function timeSeconds = convertTimeVectorToSecondsLocal( ...
    timeVector,timeUnits)

    switch timeUnits

        case 'seconds'

            timeSeconds = timeVector;

        case 'milliseconds'

            timeSeconds = timeVector./1000;

        case 'auto'

            differences = diff(timeVector);

            differences = differences( ...
                isfinite(differences) & differences>0);

            if isempty(differences)

                timeSeconds = nan(size(timeVector));
                return
            end

            medianStep = median(differences);

            if medianStep>0.5
                timeSeconds = timeVector./1000;
            else
                timeSeconds = timeVector;
            end
    end
end


% function [ellipseX,ellipseY] = covarianceEllipseLocal( ...
%     meanPoint,covarianceMatrix,ellipseScale,numberPoints)
% 
%     covarianceMatrix = ...
%         (covarianceMatrix+covarianceMatrix')./2;
% 
%     [eigenvectors,eigenvalueMatrix] = ...
%         eig(covarianceMatrix);
% 
%     eigenvalues = diag(eigenvalueMatrix);
%     eigenvalues(eigenvalues<0) = 0;
% 
%     theta = linspace(0,2*pi,numberPoints);
% 
%     unitCircle = [cos(theta);sin(theta)];
% 
%     transform = ...
%         eigenvectors * ...
%         diag(sqrt(eigenvalues)) * ...
%         ellipseScale;
% 
%     ellipsePoints = transform*unitCircle;
% 
%     ellipseX = ellipsePoints(1,:)+meanPoint(1);
%     ellipseY = ellipsePoints(2,:)+meanPoint(2);
% end


function plotPairedValuesLocal(firstValues,secondValues)

    numberAnimals = numel(firstValues);

    for animalIdx = 1:numberAnimals

        plot( ...
            [1 2], ...
            [firstValues(animalIdx),secondValues(animalIdx)], ...
            '-o', ...
            'Color',[0.65 0.65 0.65], ...
            'MarkerFaceColor','w', ...
            'LineWidth',1.2);
    end

    groupMeans = [ ...
        mean(firstValues,'omitnan'), ...
        mean(secondValues,'omitnan')];

    groupSEM = [ ...
        std(firstValues,'omitnan')./sqrt(numberAnimals), ...
        std(secondValues,'omitnan')./sqrt(numberAnimals)];

    errorbar( ...
        [1 2], ...
        groupMeans, ...
        groupSEM, ...
        'ko-', ...
        'MarkerFaceColor','k', ...
        'LineWidth',2.3);

    xlim([0.7 2.3]);
end

%% ========================================================================
% FIGURE 2: ALL INDIVIDUAL PREDICTED AND ACTUAL ENDPOINTS
%
% Each panel contains all valid reaches across all animals for one subtype.
%
% Small gray points:
%   Deceleration-aware predicted endpoints
%
% Small orange points:
%   Actual endpoints at maximum outward x position
%
% Light connecting lines:
%   Predicted-to-actual shift for the same reach
%
% Large gray point:
%   Mean predicted endpoint
%
% Large orange point:
%   Mean actual endpoint
%
% Black arrow:
%   Mean actual endpoint - mean predicted endpoint
%
% Dotted line:
%   Direction from the mean predicted endpoint to the target
%
% Black star:
%   Target
% ========================================================================

figure( ...
    'Name','Individual deceleration-aware endpoint predictions', ...
    'Color','w', ...
    'Position',[100 60 1200 950]);

layoutEndpoints = tiledlayout( ...
    2,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');


for subtypeIdx = 1:numberSubtypes

    subtype = subtypes{subtypeIdx};

    R = DecelerationEndpointResults.(subtype);

    ax = nexttile;

    hold(ax,'on');


    %% --------------------------------------------------------------------
    % Extract valid paired endpoints
    % ---------------------------------------------------------------------

    predictedEndpoints = ...
        R.reach.predictedEndpoint .* positionScale;

    actualEndpoints = ...
        R.reach.actualEndpoint .* positionScale;


    validPairs = ...
        all(isfinite(predictedEndpoints),2) & ...
        all(isfinite(actualEndpoints),2);


    predictedEndpoints = ...
        predictedEndpoints(validPairs,:);

    actualEndpoints = ...
        actualEndpoints(validPairs,:);


    numberValidReaches = ...
        size(predictedEndpoints,1);


    if numberValidReaches == 0

        title(ax,sprintf('%s: no valid endpoint data',subtype));

        xlim(ax,[1.5 4.5]);
        ylim(ax,[0.5 3.5]);

        xlabel(ax,sprintf('x position (%s)',positionUnit));
        ylabel(ax,sprintf('y position (%s)',positionUnit));

        axis(ax,'equal');
        box(ax,'off');

        continue

    end


    %% --------------------------------------------------------------------
    % Plot paired predicted-to-actual endpoint shifts
    %
    % Each line connects endpoints from the same reach.
    % ---------------------------------------------------------------------

    for reachIdx = 1:numberValidReaches

        plot( ...
            ax, ...
            [ ...
                predictedEndpoints(reachIdx,1), ...
                actualEndpoints(reachIdx,1)], ...
            [ ...
                predictedEndpoints(reachIdx,2), ...
                actualEndpoints(reachIdx,2)], ...
            '-', ...
            'Color',[1 1 1], ...
            'LineWidth',0.35);

    end


    %% --------------------------------------------------------------------
    % Plot all predicted endpoints
    % ---------------------------------------------------------------------

    predictedPointsHandle = scatter( ...
        ax, ...
        predictedEndpoints(:,1), ...
        predictedEndpoints(:,2), ...
        20, ...
        [0.35 0.35 0.35], ...
        'filled', ...
        'MarkerFaceAlpha',0.30, ...
        'MarkerEdgeAlpha',0.15);


    %% --------------------------------------------------------------------
    % Plot all actual endpoints
    % ---------------------------------------------------------------------

    actualPointsHandle = scatter( ...
        ax, ...
        actualEndpoints(:,1), ...
        actualEndpoints(:,2), ...
        20, ...
        [ 0.5781         0    0.8242], ...
        'filled', ...
        'MarkerFaceAlpha',0.30, ...
        'MarkerEdgeAlpha',0.15);


    %% --------------------------------------------------------------------
    % Calculate mean endpoints
    % ---------------------------------------------------------------------

    meanPredictedEndpoint = ...
        mean(predictedEndpoints,1,'omitnan');

    meanActualEndpoint = ...
        mean(actualEndpoints,1,'omitnan');


    %% --------------------------------------------------------------------
    % Mean correction vector
    % ---------------------------------------------------------------------

    meanCorrectionVector = ...
        meanActualEndpoint - meanPredictedEndpoint;


    correctionArrowHandle = quiver( ...
        ax, ...
        meanPredictedEndpoint(1), ...
        meanPredictedEndpoint(2), ...
        meanCorrectionVector(1), ...
        meanCorrectionVector(2), ...
        0, ...
        'Color','k', ...
        'LineWidth',3, ...
        'MaxHeadSize',1.2);


    %% --------------------------------------------------------------------
    % Plot mean endpoints
    % ---------------------------------------------------------------------

    predictedMeanHandle = scatter( ...
        ax, ...
        meanPredictedEndpoint(1), ...
        meanPredictedEndpoint(2), ...
        150, ...
        [0.35 0.35 0.35], ...
        'filled', ...
        'MarkerEdgeColor','k', ...
        'LineWidth',1.5);


    actualMeanHandle = scatter( ...
        ax, ...
        meanActualEndpoint(1), ...
        meanActualEndpoint(2), ...
        150, ...
        [0.90 0.40 0.10], ...
        'filled', ...
        'MarkerEdgeColor','k', ...
        'LineWidth',1.5);


    %% --------------------------------------------------------------------
    % Plot target
    % ---------------------------------------------------------------------

    target2DForPlot = ...
        targetPosition(analysisDimensions) .* positionScale;


    targetHandle = scatter( ...
        ax, ...
        target2DForPlot(1), ...
        target2DForPlot(2), ...
        200, ...
        'k', ...
        'pentagram', ...
        'filled');


    %% --------------------------------------------------------------------
    % Draw direction from mean predicted endpoint to target
    % ---------------------------------------------------------------------

    targetDirectionHandle = plot( ...
        ax, ...
        [ ...
            meanPredictedEndpoint(1), ...
            target2DForPlot(1)], ...
        [ ...
            meanPredictedEndpoint(2), ...
            target2DForPlot(2)], ...
        'k:', ...
        'LineWidth',1.5);


    %% --------------------------------------------------------------------
    % Formatting
    % ---------------------------------------------------------------------

    axis(ax,'equal');

    % xlim(ax,[1.5 4.5]);
    % ylim(ax,[0.5 3.5]);

    xlabel(ax,sprintf('x position (%s)',positionUnit));
    ylabel(ax,sprintf('y position (%s)',positionUnit));


    title(ax,sprintf( ...
        '%s: predicted and actual endpoints, n = %d', ...
        subtype, ...
        numberValidReaches));


    box(ax,'off');


    set( ...
        ax, ...
        'FontSize',11, ...
        'LineWidth',1.2, ...
        'TickDir','out');


    if subtypeIdx == 1

        legend( ...
            ax, ...
            [ ...
                predictedPointsHandle, ...
                actualPointsHandle, ...
                predictedMeanHandle, ...
                actualMeanHandle, ...
                correctionArrowHandle, ...
                targetDirectionHandle, ...
                targetHandle], ...
            { ...
                'Individual predicted endpoints', ...
                'Individual actual endpoints', ...
                'Mean predicted endpoint', ...
                'Mean actual endpoint', ...
                'Mean actual-minus-predicted shift', ...
                'Mean prediction-to-target direction', ...
                'Target'}, ...
            'Location','bestoutside');

    end

end


%% NEW CODE WITH EDITED FIGURE

%% ========================================================================
% FIGURE 2: INDIVIDUAL PREDICTED AND ACTUAL ENDPOINTS
% WITH TARGET-ERROR STATISTICS
%
% For each subtype:
%
%   Small gray points:
%       Deceleration-aware predicted endpoints
%
%   Small purple points:
%       Actual endpoints at maximum outward x position
%
%   Light connecting lines:
%       Predicted-to-actual shift for the same reach
%
%   Large gray point:
%       Mean predicted endpoint
%
%   Large purple point:
%       Mean actual endpoint
%
%   Black arrow:
%       Mean actual endpoint - mean predicted endpoint
%
%   Dotted line:
%       Direction from mean predicted endpoint to target
%
% STATISTICS
%
% For each reach:
%
%   squared target error =
%       (endpoint x - target x)^2 +
%       (endpoint y - target y)^2
%
% MSE is first calculated within each animal. Predicted and actual
% animal-level MSE values are then compared with a paired Wilcoxon
% signed-rank test.
%
% RMSE is also reported because it has the same units as position.
% ========================================================================


%% ------------------------------------------------------------------------
% Plot and statistical settings
% -------------------------------------------------------------------------

% These limits are applied identically to every subtype panel.
commonXLimits = [1.5 4.5];
commonYLimits = [0.5 3.5];

% Position coordinates are already multiplied by positionScale.
%
% Example:
%   positionScale = 10;
%   positionUnit = 'cm';

predictedPointColor = [0.35 0.35 0.35];
actualPointColor = [0.5781 0 0.8242];

predictedMeanColor = [0.35 0.35 0.35];
actualMeanColor = [0.5781 0 0.8242];

connectionLineColor = [0.82 0.82 0.82];

% Allocate subtype-level statistics.
SubtypeEndpointStats = struct;

predictedMSEBySubtype = nan(numberSubtypes,1);
actualMSEBySubtype = nan(numberSubtypes,1);

predictedRMSEBySubtype = nan(numberSubtypes,1);
actualRMSEBySubtype = nan(numberSubtypes,1);

pMSEBySubtype = nan(numberSubtypes,1);

numberAnimalsBySubtype = nan(numberSubtypes,1);
numberReachesBySubtype = nan(numberSubtypes,1);


%% ========================================================================
% FIRST PASS: CALCULATE REACH- AND ANIMAL-LEVEL STATISTICS
% ========================================================================

for subtypeIdx = 1:numberSubtypes

    subtype = subtypes{subtypeIdx};

    R = DecelerationEndpointResults.(subtype);


    %% Extract endpoint pairs in plotting units

    predictedEndpoints = ...
        R.reach.predictedEndpoint .* positionScale;

    actualEndpoints = ...
        R.reach.actualEndpoint .* positionScale;

    animalIDs = ...
        R.reach.animalID;


    %% Retain matched, valid endpoint pairs

    validPairs = ...
        all(isfinite(predictedEndpoints),2) & ...
        all(isfinite(actualEndpoints),2) & ...
        isfinite(animalIDs);


    predictedEndpoints = ...
        predictedEndpoints(validPairs,:);

    actualEndpoints = ...
        actualEndpoints(validPairs,:);

    animalIDs = ...
        animalIDs(validPairs);


    numberValidReaches = ...
        size(predictedEndpoints,1);

    numberReachesBySubtype(subtypeIdx) = ...
        numberValidReaches;


    if numberValidReaches == 0

        SubtypeEndpointStats.(subtype).predictedEndpoints = [];
        SubtypeEndpointStats.(subtype).actualEndpoints = [];
        SubtypeEndpointStats.(subtype).animalIDs = [];

        continue

    end


    %% Target in the same units as the plotted endpoints

    target2DForPlot = ...
        targetPosition(analysisDimensions) .* positionScale;


    %% --------------------------------------------------------------------
    % Reach-level squared endpoint errors
    %
    % squared error =
    %   x error squared + y error squared
    % ---------------------------------------------------------------------

    predictedTargetDifference = ...
        predictedEndpoints - target2DForPlot;

    actualTargetDifference = ...
        actualEndpoints - target2DForPlot;


    predictedSquaredError = ...
        sum(predictedTargetDifference.^2,2);

    actualSquaredError = ...
        sum(actualTargetDifference.^2,2);


    %% Euclidean errors, useful for interpretation

    predictedEuclideanError = ...
        sqrt(predictedSquaredError);

    actualEuclideanError = ...
        sqrt(actualSquaredError);


    %% --------------------------------------------------------------------
    % Animal-level MSE
    %
    % Each animal contributes one predicted MSE and one actual MSE.
    % ---------------------------------------------------------------------

    uniqueAnimalIDs = ...
        unique(animalIDs);

    numberAnimals = ...
        numel(uniqueAnimalIDs);


    animalPredictedMSE = ...
        nan(numberAnimals,1);

    animalActualMSE = ...
        nan(numberAnimals,1);

    animalPredictedRMSE = ...
        nan(numberAnimals,1);

    animalActualRMSE = ...
        nan(numberAnimals,1);

    animalReachCount = ...
        nan(numberAnimals,1);


    for animalIdx = 1:numberAnimals

        thisAnimalID = ...
            uniqueAnimalIDs(animalIdx);

        thisAnimalMask = ...
            animalIDs == thisAnimalID;

        animalReachCount(animalIdx) = ...
            sum(thisAnimalMask);


        animalPredictedMSE(animalIdx) = ...
            mean( ...
                predictedSquaredError(thisAnimalMask), ...
                'omitnan');


        animalActualMSE(animalIdx) = ...
            mean( ...
                actualSquaredError(thisAnimalMask), ...
                'omitnan');


        animalPredictedRMSE(animalIdx) = ...
            sqrt(animalPredictedMSE(animalIdx));


        animalActualRMSE(animalIdx) = ...
            sqrt(animalActualMSE(animalIdx));

    end


    %% Retain animals with paired values

    validAnimals = ...
        isfinite(animalPredictedMSE) & ...
        isfinite(animalActualMSE);


    animalPredictedMSE = ...
        animalPredictedMSE(validAnimals);

    animalActualMSE = ...
        animalActualMSE(validAnimals);

    animalPredictedRMSE = ...
        animalPredictedRMSE(validAnimals);

    animalActualRMSE = ...
        animalActualRMSE(validAnimals);

    validUniqueAnimalIDs = ...
        uniqueAnimalIDs(validAnimals);

    validAnimalReachCount = ...
        animalReachCount(validAnimals);


    numberValidAnimals = ...
        sum(validAnimals);

    numberAnimalsBySubtype(subtypeIdx) = ...
        numberValidAnimals;


    %% --------------------------------------------------------------------
    % Paired animal-level statistical test
    %
    % Directional alternative:
    %   actual MSE < predicted MSE
    % ---------------------------------------------------------------------

    if numberValidAnimals >= 2

        pMSE = signrank( ...
            animalActualMSE, ...
            animalPredictedMSE, ...
            'tail','left');

    else

        pMSE = NaN;

    end


    %% Group summaries based on animal-level values

    groupPredictedMSE = ...
        mean(animalPredictedMSE,'omitnan');

    groupActualMSE = ...
        mean(animalActualMSE,'omitnan');

    groupPredictedRMSE = ...
        mean(animalPredictedRMSE,'omitnan');

    groupActualRMSE = ...
        mean(animalActualRMSE,'omitnan');


    predictedMSEBySubtype(subtypeIdx) = ...
        groupPredictedMSE;

    actualMSEBySubtype(subtypeIdx) = ...
        groupActualMSE;

    predictedRMSEBySubtype(subtypeIdx) = ...
        groupPredictedRMSE;

    actualRMSEBySubtype(subtypeIdx) = ...
        groupActualRMSE;

    pMSEBySubtype(subtypeIdx) = ...
        pMSE;


    %% --------------------------------------------------------------------
    % Store all values
    % ---------------------------------------------------------------------

    SubtypeEndpointStats.(subtype).predictedEndpoints = ...
        predictedEndpoints;

    SubtypeEndpointStats.(subtype).actualEndpoints = ...
        actualEndpoints;

    SubtypeEndpointStats.(subtype).animalIDs = ...
        animalIDs;

    SubtypeEndpointStats.(subtype).target = ...
        target2DForPlot;


    SubtypeEndpointStats.(subtype).reach.predictedSquaredError = ...
        predictedSquaredError;

    SubtypeEndpointStats.(subtype).reach.actualSquaredError = ...
        actualSquaredError;

    SubtypeEndpointStats.(subtype).reach.predictedEuclideanError = ...
        predictedEuclideanError;

    SubtypeEndpointStats.(subtype).reach.actualEuclideanError = ...
        actualEuclideanError;


    SubtypeEndpointStats.(subtype).animal.ID = ...
        validUniqueAnimalIDs;

    SubtypeEndpointStats.(subtype).animal.reachCount = ...
        validAnimalReachCount;

    SubtypeEndpointStats.(subtype).animal.predictedMSE = ...
        animalPredictedMSE;

    SubtypeEndpointStats.(subtype).animal.actualMSE = ...
        animalActualMSE;

    SubtypeEndpointStats.(subtype).animal.predictedRMSE = ...
        animalPredictedRMSE;

    SubtypeEndpointStats.(subtype).animal.actualRMSE = ...
        animalActualRMSE;


    SubtypeEndpointStats.(subtype).group.predictedMSE = ...
        groupPredictedMSE;

    SubtypeEndpointStats.(subtype).group.actualMSE = ...
        groupActualMSE;

    SubtypeEndpointStats.(subtype).group.predictedRMSE = ...
        groupPredictedRMSE;

    SubtypeEndpointStats.(subtype).group.actualRMSE = ...
        groupActualRMSE;

    SubtypeEndpointStats.(subtype).statistics.pMSE = ...
        pMSE;


    %% Print subtype results

    fprintf('\n%s endpoint-target error\n',subtype);

    fprintf('Animals: %d\n',numberValidAnimals);
    fprintf('Reaches: %d\n',numberValidReaches);

    fprintf('Predicted MSE: %.5f %s^2\n', ...
        groupPredictedMSE,positionUnit);

    fprintf('Actual MSE: %.5f %s^2\n', ...
        groupActualMSE,positionUnit);

    fprintf('Predicted RMSE: %.5f %s\n', ...
        groupPredictedRMSE,positionUnit);

    fprintf('Actual RMSE: %.5f %s\n', ...
        groupActualRMSE,positionUnit);

    fprintf('Paired animal-level MSE p = %.6g\n', ...
        pMSE);

end


%% ========================================================================
% OPTIONAL FDR CORRECTION ACROSS T0-T3
% ========================================================================

qMSEBySubtype = ...
    benjaminiHochbergLocal(pMSEBySubtype);


for subtypeIdx = 1:numberSubtypes

    subtype = subtypes{subtypeIdx};

    SubtypeEndpointStats.(subtype).statistics.qMSE = ...
        qMSEBySubtype(subtypeIdx);

end


%% ========================================================================
% FIGURE: ALL INDIVIDUAL ENDPOINTS
% ========================================================================

figureEndpointScatter = figure( ...
    'Name','Individual deceleration-aware endpoint predictions', ...
    'Color','w', ...
    'Position',[100 60 1250 950]);


layoutEndpoints = tiledlayout( ...
    figureEndpointScatter, ...
    2,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');


% Store handles for one shared legend.
sharedLegendHandles = gobjects(7,1);


for subtypeIdx = 1:numberSubtypes

    subtype = subtypes{subtypeIdx};

    S = SubtypeEndpointStats.(subtype);


    ax = nexttile(layoutEndpoints,subtypeIdx);

    hold(ax,'on');


    %% --------------------------------------------------------------------
    % No-data panel
    % ---------------------------------------------------------------------

    if isempty(S.predictedEndpoints)

        title(ax,sprintf('%s: no valid endpoint data',subtype));

        xlim(ax,commonXLimits);
        ylim(ax,commonYLimits);

        xlabel(ax,sprintf('x position (%s)',positionUnit));
        ylabel(ax,sprintf('y position (%s)',positionUnit));

        axis(ax,'equal');
        pbaspect(ax,[1 1 1]);

        box(ax,'off');

        continue

    end


    predictedEndpoints = ...
        S.predictedEndpoints;

    actualEndpoints = ...
        S.actualEndpoints;

    target2DForPlot = ...
        S.target;

    numberValidReaches = ...
        size(predictedEndpoints,1);


    %% --------------------------------------------------------------------
    % Plot paired predicted-to-actual shifts
    % ---------------------------------------------------------------------

    for reachIdx = 1:numberValidReaches

        plot( ...
            ax, ...
            [ ...
                predictedEndpoints(reachIdx,1), ...
                actualEndpoints(reachIdx,1)], ...
            [ ...
                predictedEndpoints(reachIdx,2), ...
                actualEndpoints(reachIdx,2)], ...
            '-', ...
            'Color',[1 1 1], ...
            'LineWidth',0.35, ...
            'HandleVisibility','off');

    end


    %% --------------------------------------------------------------------
    % Individual predicted endpoints
    % ---------------------------------------------------------------------

    predictedPointsHandle = scatter( ...
        ax, ...
        predictedEndpoints(:,1), ...
        predictedEndpoints(:,2), ...
        20, ...
        predictedPointColor, ...
        'filled', ...
        'MarkerFaceAlpha',0.30, ...
        'MarkerEdgeAlpha',0.15);


    %% --------------------------------------------------------------------
    % Individual actual endpoints
    % ---------------------------------------------------------------------

    actualPointsHandle = scatter( ...
        ax, ...
        actualEndpoints(:,1), ...
        actualEndpoints(:,2), ...
        20, ...
        actualPointColor, ...
        'filled', ...
        'MarkerFaceAlpha',0.30, ...
        'MarkerEdgeAlpha',0.15);


    %% --------------------------------------------------------------------
    % Mean endpoints
    % ---------------------------------------------------------------------

    meanPredictedEndpoint = ...
        mean(predictedEndpoints,1,'omitnan');

    meanActualEndpoint = ...
        mean(actualEndpoints,1,'omitnan');


    meanCorrectionVector = ...
        meanActualEndpoint - meanPredictedEndpoint;


    correctionArrowHandle = quiver( ...
        ax, ...
        meanPredictedEndpoint(1), ...
        meanPredictedEndpoint(2), ...
        meanCorrectionVector(1), ...
        meanCorrectionVector(2), ...
        0, ...
        'Color','k', ...
        'LineWidth',3, ...
        'MaxHeadSize',1.2);


    predictedMeanHandle = scatter( ...
        ax, ...
        meanPredictedEndpoint(1), ...
        meanPredictedEndpoint(2), ...
        150, ...
        predictedMeanColor, ...
        'filled', ...
        'MarkerEdgeColor','k', ...
        'LineWidth',1.5);


    actualMeanHandle = scatter( ...
        ax, ...
        meanActualEndpoint(1), ...
        meanActualEndpoint(2), ...
        150, ...
        actualMeanColor, ...
        'filled', ...
        'MarkerEdgeColor','k', ...
        'LineWidth',1.5);


    %% --------------------------------------------------------------------
    % Target and target direction
    % ---------------------------------------------------------------------

    targetHandle = scatter( ...
        ax, ...
        target2DForPlot(1), ...
        target2DForPlot(2), ...
        200, ...
        'k', ...
        'pentagram', ...
        'filled');


    targetDirectionHandle = plot( ...
        ax, ...
        [ ...
            meanPredictedEndpoint(1), ...
            target2DForPlot(1)], ...
        [ ...
            meanPredictedEndpoint(2), ...
            target2DForPlot(2)], ...
        'k:', ...
        'LineWidth',1.5);


    %% --------------------------------------------------------------------
    % Formatting: same limits and same panel geometry
    % ---------------------------------------------------------------------

    % xlim(ax,commonXLimits);
    % ylim(ax,commonYLimits);
    xlim([0 4]);
    ylim([0 4]);

    axis(ax,'equal');

    % Forces each data panel to have the same square plotting-box shape.
    pbaspect(ax,[1 1 1]);


    xlabel(ax,sprintf('x position (%s)',positionUnit));
    ylabel(ax,sprintf('y position (%s)',positionUnit));


    title(ax,{ ...
        sprintf('%s: n = %d reaches, N = %d animals', ...
            subtype, ...
            numberValidReaches, ...
            numberAnimalsBySubtype(subtypeIdx)), ...
        sprintf('RMSE: predicted %.3f, actual %.3f %s; p = %.3g', ...
            predictedRMSEBySubtype(subtypeIdx), ...
            actualRMSEBySubtype(subtypeIdx), ...
            positionUnit, ...
            pMSEBySubtype(subtypeIdx))});


    box(ax,'off');


    set( ...
        ax, ...
        'FontSize',11, ...
        'LineWidth',1.2, ...
        'TickDir','out');


    %% Save handles from the first populated panel for a shared legend

    if all(~isgraphics(sharedLegendHandles))

        sharedLegendHandles = [ ...
            predictedPointsHandle; ...
            actualPointsHandle; ...
            predictedMeanHandle; ...
            actualMeanHandle; ...
            correctionArrowHandle; ...
            targetDirectionHandle; ...
            targetHandle];

    end

end


%% Shared legend
%
% Putting the legend in the tiled layout rather than inside one axis keeps
% all four endpoint panels the same size.

if all(isgraphics(sharedLegendHandles))

    endpointLegend = legend( ...
        sharedLegendHandles, ...
        { ...
            'Individual predicted endpoints', ...
            'Individual actual endpoints', ...
            'Mean predicted endpoint', ...
            'Mean actual endpoint', ...
            'Mean actual-minus-predicted shift', ...
            'Mean prediction-to-target direction', ...
            'Target'}, ...
        'Orientation','horizontal', ...
        'NumColumns',4);

    endpointLegend.Layout.Tile = 'south';

end


%% ========================================================================
% FIGURE: PAIRED ANIMAL-LEVEL MSE
%
% This is the inferential companion to the endpoint scatter plot.
% ========================================================================

figureMSE = figure( ...
    'Name','Predicted versus actual endpoint MSE', ...
    'Color','w', ...
    'Position',[150 100 1100 850]);


layoutMSE = tiledlayout( ...
    figureMSE, ...
    2,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');


for subtypeIdx = 1:numberSubtypes

    subtype = subtypes{subtypeIdx};

    S = SubtypeEndpointStats.(subtype);

    ax = nexttile(layoutMSE,subtypeIdx);

    hold(ax,'on');


    if ~isfield(S,'animal') || ...
            isempty(S.animal.predictedMSE)

        title(ax,sprintf('%s: insufficient animal data',subtype));

        box(ax,'off');

        continue

    end


    predictedAnimalMSE = ...
        S.animal.predictedMSE;

    actualAnimalMSE = ...
        S.animal.actualMSE;

    numberAnimals = ...
        numel(predictedAnimalMSE);


    %% Plot every animal as a paired line

    for animalIdx = 1:numberAnimals

        plot( ...
            ax, ...
            [1 2], ...
            [ ...
                predictedAnimalMSE(animalIdx), ...
                actualAnimalMSE(animalIdx)], ...
            '-o', ...
            'Color',[0.65 0.65 0.65], ...
            'MarkerFaceColor','w', ...
            'LineWidth',1.25, ...
            'MarkerSize',6);

    end


    %% Group mean and SEM

    groupMeans = [ ...
        mean(predictedAnimalMSE,'omitnan'), ...
        mean(actualAnimalMSE,'omitnan')];


    groupSEM = [ ...
        std(predictedAnimalMSE,'omitnan') ./ ...
            sqrt(numberAnimals), ...
        std(actualAnimalMSE,'omitnan') ./ ...
            sqrt(numberAnimals)];


    errorbar( ...
        ax, ...
        [1 2], ...
        groupMeans, ...
        groupSEM, ...
        'ko-', ...
        'MarkerFaceColor','k', ...
        'LineWidth',2.4, ...
        'MarkerSize',8);


    xlim(ax,[0.7 2.3]);

    xticks(ax,[1 2]);

    xticklabels(ax,{ ...
        'Predicted', ...
        'Actual'});


    ylabel(ax,sprintf('Mean squared endpoint error (%s^2)', ...
        positionUnit));


    title(ax,sprintf( ...
        '%s: p = %.3g, FDR q = %.3g', ...
        subtype, ...
        pMSEBySubtype(subtypeIdx), ...
        qMSEBySubtype(subtypeIdx)));


    box(ax,'off');


    set( ...
        ax, ...
        'FontSize',11, ...
        'LineWidth',1.2, ...
        'TickDir','out');

end

%% ========================================================================
% ALL-CSM DECELERATION-AWARE ENDPOINT ANALYSIS
%
% Combines valid reaches from:
%
%   T0
%   T1
%   T2
%   T3
%
% into a single "All CSMs" category.
%
% Importantly:
%
%   - Reach subtype is ignored for the primary combined analysis.
%   - Animal identity is retained.
%   - Statistics are performed on animal-level means, not pooled reaches.
%
% REQUIRED INPUT:
%
%   DecelerationEndpointResults.T0
%   DecelerationEndpointResults.T1
%   DecelerationEndpointResults.T2
%   DecelerationEndpointResults.T3
%
% Each subtype must contain:
%
%   R.reach.predictedEndpoint
%   R.reach.actualEndpoint
%   R.reach.animalID
% ========================================================================


%% ------------------------------------------------------------------------
% Settings
% -------------------------------------------------------------------------

combinedSubtypes = {'T0','T1','T2','T3'};

% Same plotting units used in the preceding analysis.
%
% For centimeters:
%   positionScale = 10;
%   positionUnit  = 'cm';
%
% For decimeters:
%   positionScale = 1;
%   positionUnit  = 'dm';

positionScale = 10;
positionUnit = 'cm';

% Apply the same limits used in the subtype plots.
commonXLimits = [1.5 4.5];
commonYLimits = [0.5 3.5];

% Colors
predictedPointColor = [0.35 0.35 0.35];
actualPointColor = [0.5781 0 0.8242];

connectionLineColor = [0.82 0.82 0.82];

% Minimum number of valid pooled reaches required for an animal to
% contribute to the animal-level statistical analysis.
minimumCombinedReachesPerAnimal = 5;


%% ========================================================================
% COMBINE ALL VALID REACHES
% ========================================================================

allPredictedEndpoints = [];
allActualEndpoints = [];

allAnimalIDs = [];
allSubtypeIDs = [];
allSubtypeLabels = {};

allOriginalReachIndices = [];


for subtypeIdx = 1:numel(combinedSubtypes)

    subtype = combinedSubtypes{subtypeIdx};

    if ~isfield(DecelerationEndpointResults,subtype)

        warning('%s is absent from DecelerationEndpointResults.',subtype);
        continue

    end

    R = DecelerationEndpointResults.(subtype);


    %% Confirm required fields exist

    requiredFieldsPresent = ...
        isfield(R,'reach') && ...
        isfield(R.reach,'predictedEndpoint') && ...
        isfield(R.reach,'actualEndpoint') && ...
        isfield(R.reach,'animalID');


    if ~requiredFieldsPresent

        warning('%s does not contain the required reach-level fields.', ...
            subtype);

        continue

    end


    %% Extract subtype data in original units

    predictedEndpoints = ...
        R.reach.predictedEndpoint;

    actualEndpoints = ...
        R.reach.actualEndpoint;

    animalIDs = ...
        R.reach.animalID;


    %% Original reach index is optional

    if isfield(R.reach,'originalReachIndex')

        originalReachIndices = ...
            R.reach.originalReachIndex;

    else

        originalReachIndices = ...
            (1:size(predictedEndpoints,1))';

    end


    %% Retain matched, valid endpoint pairs

    validPairs = ...
        all(isfinite(predictedEndpoints),2) & ...
        all(isfinite(actualEndpoints),2) & ...
        isfinite(animalIDs);


    predictedEndpoints = ...
        predictedEndpoints(validPairs,:);

    actualEndpoints = ...
        actualEndpoints(validPairs,:);

    animalIDs = ...
        animalIDs(validPairs);

    originalReachIndices = ...
        originalReachIndices(validPairs);


    numberSubtypeReaches = ...
        size(predictedEndpoints,1);


    %% Add subtype reaches to the combined arrays

    allPredictedEndpoints = [ ...
        allPredictedEndpoints; ...
        predictedEndpoints];

    allActualEndpoints = [ ...
        allActualEndpoints; ...
        actualEndpoints];

    allAnimalIDs = [ ...
        allAnimalIDs; ...
        animalIDs];

    allSubtypeIDs = [ ...
        allSubtypeIDs; ...
        repmat(subtypeIdx,numberSubtypeReaches,1)];

    allSubtypeLabels = [ ...
        allSubtypeLabels; ...
        repmat({subtype},numberSubtypeReaches,1)];

    allOriginalReachIndices = [ ...
        allOriginalReachIndices; ...
        originalReachIndices];


    fprintf('%s: added %d valid reaches\n', ...
        subtype,numberSubtypeReaches);

end


numberCombinedReaches = ...
    size(allPredictedEndpoints,1);


if numberCombinedReaches == 0

    error('No valid T0-T3 reaches were available for the combined analysis.');

end


fprintf('\nTotal combined CSM reaches: %d\n', ...
    numberCombinedReaches);


%% ========================================================================
% CALCULATE REACH-LEVEL TARGET ERRORS
%
% Keep calculations in original position units.
% ========================================================================

target2D = ...
    targetPosition(analysisDimensions);


predictedTargetDifference = ...
    allPredictedEndpoints - target2D;

actualTargetDifference = ...
    allActualEndpoints - target2D;


%% Squared Euclidean target error

predictedSquaredError = ...
    sum(predictedTargetDifference.^2,2);

actualSquaredError = ...
    sum(actualTargetDifference.^2,2);


%% Euclidean target error

predictedEuclideanError = ...
    sqrt(predictedSquaredError);

actualEuclideanError = ...
    sqrt(actualSquaredError);


%% Improvement
%
% Positive:
%   actual endpoint is closer to target than predicted endpoint.

endpointErrorImprovement = ...
    predictedEuclideanError - actualEuclideanError;


%% Squared-error improvement
%
% Positive:
%   actual squared error is lower.

squaredErrorImprovement = ...
    predictedSquaredError - actualSquaredError;


%% Fraction of pooled reaches that improve

combinedFractionImprovedReachLevel = ...
    mean(endpointErrorImprovement > 0,'omitnan');


%% ========================================================================
% ANIMAL-LEVEL COMBINED ANALYSIS
% ========================================================================

uniqueAnimalIDs = ...
    unique(allAnimalIDs);

numberCandidateAnimals = ...
    numel(uniqueAnimalIDs);


animalCombinedReachCount = ...
    zeros(numberCandidateAnimals,1);

animalPredictedMSE = ...
    nan(numberCandidateAnimals,1);

animalActualMSE = ...
    nan(numberCandidateAnimals,1);

animalPredictedRMSE = ...
    nan(numberCandidateAnimals,1);

animalActualRMSE = ...
    nan(numberCandidateAnimals,1);

animalMeanPredictedError = ...
    nan(numberCandidateAnimals,1);

animalMeanActualError = ...
    nan(numberCandidateAnimals,1);

animalMeanImprovement = ...
    nan(numberCandidateAnimals,1);

animalFractionImproved = ...
    nan(numberCandidateAnimals,1);

% Count T0-T3 reaches separately for each animal.
animalSubtypeReachCounts = ...
    zeros(numberCandidateAnimals,numel(combinedSubtypes));


for animalIdx = 1:numberCandidateAnimals

    thisAnimalID = ...
        uniqueAnimalIDs(animalIdx);

    animalMask = ...
        allAnimalIDs == thisAnimalID;


    animalCombinedReachCount(animalIdx) = ...
        sum(animalMask);


    %% Count subtype contributions

    for subtypeIdx = 1:numel(combinedSubtypes)

        animalSubtypeReachCounts(animalIdx,subtypeIdx) = ...
            sum( ...
                animalMask & ...
                allSubtypeIDs == subtypeIdx);

    end


    %% Only summarize animals with enough combined reaches

    if animalCombinedReachCount(animalIdx) < ...
            minimumCombinedReachesPerAnimal

        continue

    end


    %% MSE: mean squared two-dimensional distance from target

    animalPredictedMSE(animalIdx) = ...
        mean( ...
            predictedSquaredError(animalMask), ...
            'omitnan');

    animalActualMSE(animalIdx) = ...
        mean( ...
            actualSquaredError(animalMask), ...
            'omitnan');


    %% RMSE

    animalPredictedRMSE(animalIdx) = ...
        sqrt(animalPredictedMSE(animalIdx));

    animalActualRMSE(animalIdx) = ...
        sqrt(animalActualMSE(animalIdx));


    %% Mean Euclidean target error

    animalMeanPredictedError(animalIdx) = ...
        mean( ...
            predictedEuclideanError(animalMask), ...
            'omitnan');

    animalMeanActualError(animalIdx) = ...
        mean( ...
            actualEuclideanError(animalMask), ...
            'omitnan');


    %% Mean improvement

    animalMeanImprovement(animalIdx) = ...
        mean( ...
            endpointErrorImprovement(animalMask), ...
            'omitnan');


    %% Fraction of this animal's reaches that improved

    animalFractionImproved(animalIdx) = ...
        mean( ...
            endpointErrorImprovement(animalMask) > 0, ...
            'omitnan');

end


%% Retain animals with valid paired values

validAnimalsMSE = ...
    isfinite(animalPredictedMSE) & ...
    isfinite(animalActualMSE);

validAnimalsError = ...
    isfinite(animalMeanPredictedError) & ...
    isfinite(animalMeanActualError);

validAnimalsImprovement = ...
    isfinite(animalMeanImprovement);

validAnimalsFraction = ...
    isfinite(animalFractionImproved);


numberValidAnimals = ...
    sum(validAnimalsMSE);


if numberValidAnimals < 2

    warning(['Fewer than two animals met the combined reach threshold. ' ...
        'Statistical p-values will be NaN.']);

end


%% ========================================================================
% PAIRED ANIMAL-LEVEL STATISTICS
% ========================================================================

%% Primary test: actual MSE lower than predicted MSE

if sum(validAnimalsMSE) >= 2

    pCombinedMSE = signrank( ...
        animalActualMSE(validAnimalsMSE), ...
        animalPredictedMSE(validAnimalsMSE), ...
        'tail','left');

else

    pCombinedMSE = NaN;

end


%% Mean Euclidean error: actual lower than predicted

if sum(validAnimalsError) >= 2

    pCombinedEuclideanError = signrank( ...
        animalMeanActualError(validAnimalsError), ...
        animalMeanPredictedError(validAnimalsError), ...
        'tail','left');

else

    pCombinedEuclideanError = NaN;

end


%% Endpoint improvement greater than zero

if sum(validAnimalsImprovement) >= 2

    pCombinedImprovement = signrank( ...
        animalMeanImprovement(validAnimalsImprovement), ...
        0, ...
        'tail','right');

else

    pCombinedImprovement = NaN;

end


%% Fraction of improved reaches greater than 0.5

if sum(validAnimalsFraction) >= 2

    pCombinedFractionImproved = signrank( ...
        animalFractionImproved(validAnimalsFraction), ...
        0.5, ...
        'tail','right');

else

    pCombinedFractionImproved = NaN;

end


%% ========================================================================
% GROUP SUMMARIES
%
% Group means give each animal equal weight.
% ========================================================================

groupPredictedMSE = ...
    mean(animalPredictedMSE(validAnimalsMSE),'omitnan');

groupActualMSE = ...
    mean(animalActualMSE(validAnimalsMSE),'omitnan');


groupPredictedRMSE = ...
    mean(animalPredictedRMSE(validAnimalsMSE),'omitnan');

groupActualRMSE = ...
    mean(animalActualRMSE(validAnimalsMSE),'omitnan');


groupPredictedEuclideanError = ...
    mean( ...
        animalMeanPredictedError(validAnimalsError), ...
        'omitnan');

groupActualEuclideanError = ...
    mean( ...
        animalMeanActualError(validAnimalsError), ...
        'omitnan');


groupEndpointImprovement = ...
    mean( ...
        animalMeanImprovement(validAnimalsImprovement), ...
        'omitnan');


groupFractionImproved = ...
    mean( ...
        animalFractionImproved(validAnimalsFraction), ...
        'omitnan');


%% ========================================================================
% PRINT COMBINED RESULTS
% ========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('ALL CSMs COMBINED\n');
fprintf('============================================================\n');

fprintf('Total valid reaches: %d\n', ...
    numberCombinedReaches);

fprintf('Animals passing threshold: %d\n', ...
    numberValidAnimals);


fprintf('\nReach contribution by subtype\n');

for subtypeIdx = 1:numel(combinedSubtypes)

    fprintf('%s: %d reaches\n', ...
        combinedSubtypes{subtypeIdx}, ...
        sum(allSubtypeIDs == subtypeIdx));

end


fprintf('\nTarget MSE\n');

fprintf('Predicted MSE: %.6f %s^2\n', ...
    groupPredictedMSE .* positionScale.^2, ...
    positionUnit);

fprintf('Actual MSE: %.6f %s^2\n', ...
    groupActualMSE .* positionScale.^2, ...
    positionUnit);

fprintf('Paired one-sided MSE p = %.6g\n', ...
    pCombinedMSE);


fprintf('\nTarget RMSE\n');

fprintf('Predicted RMSE: %.6f %s\n', ...
    groupPredictedRMSE .* positionScale, ...
    positionUnit);

fprintf('Actual RMSE: %.6f %s\n', ...
    groupActualRMSE .* positionScale, ...
    positionUnit);


fprintf('\nMean Euclidean target error\n');

fprintf('Predicted error: %.6f %s\n', ...
    groupPredictedEuclideanError .* positionScale, ...
    positionUnit);

fprintf('Actual error: %.6f %s\n', ...
    groupActualEuclideanError .* positionScale, ...
    positionUnit);

fprintf('Paired one-sided error p = %.6g\n', ...
    pCombinedEuclideanError);


fprintf('\nEndpoint improvement\n');

fprintf('Mean improvement: %.6f %s\n', ...
    groupEndpointImprovement .* positionScale, ...
    positionUnit);

fprintf('Improvement greater than zero p = %.6g\n', ...
    pCombinedImprovement);


fprintf('\nFraction improved\n');

fprintf('Mean animal fraction improved: %.4f\n', ...
    groupFractionImproved);

fprintf('Pooled reach fraction improved: %.4f\n', ...
    combinedFractionImprovedReachLevel);

fprintf('Animal fraction greater than 0.5 p = %.6g\n', ...
    pCombinedFractionImproved);


%% ========================================================================
% FIGURE 1: ALL INDIVIDUAL ENDPOINTS COMBINED
% ========================================================================

predictedEndpointsForPlot = ...
    allPredictedEndpoints .* positionScale;

actualEndpointsForPlot = ...
    allActualEndpoints .* positionScale;

targetForPlot = ...
    target2D .* positionScale;


figureCombinedScatter = figure( ...
    'Name','All CSM endpoint predictions combined', ...
    'Color','w', ...
    'Position',[100 80 850 800]);


axCombined = axes(figureCombinedScatter);

hold(axCombined,'on');


%% Plot paired predicted-to-actual shifts

for reachIdx = 1:numberCombinedReaches

    plot( ...
        axCombined, ...
        [ ...
            predictedEndpointsForPlot(reachIdx,1), ...
            actualEndpointsForPlot(reachIdx,1)], ...
        [ ...
            predictedEndpointsForPlot(reachIdx,2), ...
            actualEndpointsForPlot(reachIdx,2)], ...
        '-', ...
        'Color',[1 1 1], ...
        'LineWidth',0.30, ...
        'HandleVisibility','off');

end


%% Individual predicted endpoints

predictedPointsHandle = scatter( ...
    axCombined, ...
    predictedEndpointsForPlot(:,1), ...
    predictedEndpointsForPlot(:,2), ...
    20, ...
    predictedPointColor, ...
    'filled', ...
    'MarkerFaceAlpha',0.25, ...
    'MarkerEdgeAlpha',0.10);


%% Individual actual endpoints

actualPointsHandle = scatter( ...
    axCombined, ...
    actualEndpointsForPlot(:,1), ...
    actualEndpointsForPlot(:,2), ...
    20, ...
    actualPointColor, ...
    'filled', ...
    'MarkerFaceAlpha',0.25, ...
    'MarkerEdgeAlpha',0.10);


%% Mean endpoints across all reaches for descriptive visualization

meanPredictedEndpoint = ...
    mean(predictedEndpointsForPlot,1,'omitnan');

meanActualEndpoint = ...
    mean(actualEndpointsForPlot,1,'omitnan');

meanCorrectionVector = ...
    meanActualEndpoint - meanPredictedEndpoint;


correctionArrowHandle = quiver( ...
    axCombined, ...
    meanPredictedEndpoint(1), ...
    meanPredictedEndpoint(2), ...
    meanCorrectionVector(1), ...
    meanCorrectionVector(2), ...
    0, ...
    'Color','k', ...
    'LineWidth',3, ...
    'MaxHeadSize',1.2);


predictedMeanHandle = scatter( ...
    axCombined, ...
    meanPredictedEndpoint(1), ...
    meanPredictedEndpoint(2), ...
    160, ...
    predictedPointColor, ...
    'filled', ...
    'MarkerEdgeColor','k', ...
    'LineWidth',1.5);


actualMeanHandle = scatter( ...
    axCombined, ...
    meanActualEndpoint(1), ...
    meanActualEndpoint(2), ...
    160, ...
    actualPointColor, ...
    'filled', ...
    'MarkerEdgeColor','k', ...
    'LineWidth',1.5);


%% Target

targetHandle = scatter( ...
    axCombined, ...
    targetForPlot(1), ...
    targetForPlot(2), ...
    220, ...
    'k', ...
    'pentagram', ...
    'filled');


%% Predicted mean to target

targetDirectionHandle = plot( ...
    axCombined, ...
    [meanPredictedEndpoint(1),targetForPlot(1)], ...
    [meanPredictedEndpoint(2),targetForPlot(2)], ...
    'k:', ...
    'LineWidth',1.5);


%% Formatting

% xlim(axCombined,commonXLimits);
% ylim(axCombined,commonYLimits);

xlim([0 4]);
ylim([0 4]);

axis(axCombined,'equal');
pbaspect(axCombined,[1 1 1]);

xlabel(axCombined,sprintf('x position (%s)',positionUnit));
ylabel(axCombined,sprintf('y position (%s)',positionUnit));


title(axCombined,{ ...
    sprintf('All CSMs combined: n = %d reaches, N = %d animals', ...
        numberCombinedReaches,numberValidAnimals), ...
    sprintf('RMSE: predicted %.3f, actual %.3f %s; p = %.3g', ...
        groupPredictedRMSE .* positionScale, ...
        groupActualRMSE .* positionScale, ...
        positionUnit, ...
        pCombinedMSE)});


legend( ...
    axCombined, ...
    [ ...
        predictedPointsHandle, ...
        actualPointsHandle, ...
        predictedMeanHandle, ...
        actualMeanHandle, ...
        correctionArrowHandle, ...
        targetDirectionHandle, ...
        targetHandle], ...
    { ...
        'Individual predicted endpoints', ...
        'Individual actual endpoints', ...
        'Mean predicted endpoint', ...
        'Mean actual endpoint', ...
        'Mean actual-minus-predicted shift', ...
        'Mean prediction-to-target direction', ...
        'Target'}, ...
    'Location','bestoutside');


box(axCombined,'off');

set( ...
    axCombined, ...
    'FontSize',11, ...
    'LineWidth',1.2, ...
    'TickDir','out');


%% ========================================================================
% FIGURE 2: PAIRED ANIMAL-LEVEL MSE
% ========================================================================

figureCombinedMSE = figure( ...
    'Name','All CSM predicted versus actual MSE', ...
    'Color','w', ...
    'Position',[180 120 600 600]);


axMSE = axes(figureCombinedMSE);

hold(axMSE,'on');


predictedMSEForPlot = ...
    animalPredictedMSE(validAnimalsMSE) .* positionScale.^2;

actualMSEForPlot = ...
    animalActualMSE(validAnimalsMSE) .* positionScale.^2;


for animalIdx = 1:numberValidAnimals

    plot( ...
        axMSE, ...
        [1 2], ...
        [ ...
            predictedMSEForPlot(animalIdx), ...
            actualMSEForPlot(animalIdx)], ...
        '-o', ...
        'Color',[0.65 0.65 0.65], ...
        'MarkerFaceColor','w', ...
        'LineWidth',1.25, ...
        'MarkerSize',7);

end


%% Group mean and SEM

mseGroupMeans = [ ...
    mean(predictedMSEForPlot,'omitnan'), ...
    mean(actualMSEForPlot,'omitnan')];


mseGroupSEM = [ ...
    std(predictedMSEForPlot,'omitnan') ./ ...
        sqrt(numberValidAnimals), ...
    std(actualMSEForPlot,'omitnan') ./ ...
        sqrt(numberValidAnimals)];


errorbar( ...
    axMSE, ...
    [1 2], ...
    mseGroupMeans, ...
    mseGroupSEM, ...
    'ko-', ...
    'MarkerFaceColor','k', ...
    'LineWidth',2.5, ...
    'MarkerSize',9);


xlim(axMSE,[0.7 2.3]);

xticks(axMSE,[1 2]);
xticklabels(axMSE,{'Predicted','Actual'});

ylabel(axMSE,sprintf( ...
    'Mean squared endpoint error (%s^2)', ...
    positionUnit));


title(axMSE,sprintf( ...
    'All CSMs combined: N = %d, p = %.4g', ...
    numberValidAnimals, ...
    pCombinedMSE));


box(axMSE,'off');

set( ...
    axMSE, ...
    'FontSize',11, ...
    'LineWidth',1.2, ...
    'TickDir','out');


%% ========================================================================
% FIGURE 3: PAIRED ANIMAL-LEVEL RMSE
%
% Same comparison as MSE, but displayed in position units.
% The statistical test remains the MSE test above.
% ========================================================================

figureCombinedRMSE = figure( ...
    'Name','All CSM predicted versus actual RMSE', ...
    'Color','w', ...
    'Position',[220 140 600 600]);


axRMSE = axes(figureCombinedRMSE);

hold(axRMSE,'on');


predictedRMSEForPlot = ...
    animalPredictedRMSE(validAnimalsMSE) .* positionScale;

actualRMSEForPlot = ...
    animalActualRMSE(validAnimalsMSE) .* positionScale;


for animalIdx = 1:numberValidAnimals

    plot( ...
        axRMSE, ...
        [1 2], ...
        [ ...
            predictedRMSEForPlot(animalIdx), ...
            actualRMSEForPlot(animalIdx)], ...
        '-o', ...
        'Color',[0.65 0.65 0.65], ...
        'MarkerFaceColor','w', ...
        'LineWidth',1.25, ...
        'MarkerSize',7);

end


rmseGroupMeans = [ ...
    mean(predictedRMSEForPlot,'omitnan'), ...
    mean(actualRMSEForPlot,'omitnan')];


rmseGroupSEM = [ ...
    std(predictedRMSEForPlot,'omitnan') ./ ...
        sqrt(numberValidAnimals), ...
    std(actualRMSEForPlot,'omitnan') ./ ...
        sqrt(numberValidAnimals)];


errorbar( ...
    axRMSE, ...
    [1 2], ...
    rmseGroupMeans, ...
    rmseGroupSEM, ...
    'ko-', ...
    'MarkerFaceColor','k', ...
    'LineWidth',2.5, ...
    'MarkerSize',9);


xlim(axRMSE,[0.7 2.3]);

xticks(axRMSE,[1 2]);
xticklabels(axRMSE,{'Predicted','Actual'});

ylabel(axRMSE,sprintf('Root mean squared endpoint error (%s)', ...
    positionUnit));


title(axRMSE,sprintf( ...
    'All CSMs combined: N = %d', ...
    numberValidAnimals));


box(axRMSE,'off');

set( ...
    axRMSE, ...
    'FontSize',11, ...
    'LineWidth',1.2, ...
    'TickDir','out');


%% ========================================================================
% ANIMAL CONTRIBUTION TABLE
% ========================================================================

includedAnimalIDs = ...
    uniqueAnimalIDs;

if exist('animalNames','var') && ...
        numel(animalNames) >= max(includedAnimalIDs)

    Animal = ...
        animalNames(includedAnimalIDs);

else

    Animal = ...
        compose('Animal_%d',includedAnimalIDs);

end


T0Reaches = animalSubtypeReachCounts(:,1);
T1Reaches = animalSubtypeReachCounts(:,2);
T2Reaches = animalSubtypeReachCounts(:,3);
T3Reaches = animalSubtypeReachCounts(:,4);

TotalReaches = animalCombinedReachCount;

IncludedInStatistics = validAnimalsMSE;

PredictedMSE = ...
    animalPredictedMSE .* positionScale.^2;

ActualMSE = ...
    animalActualMSE .* positionScale.^2;

PredictedRMSE = ...
    animalPredictedRMSE .* positionScale;

ActualRMSE = ...
    animalActualRMSE .* positionScale;

MeanImprovement = ...
    animalMeanImprovement .* positionScale;

FractionImproved = ...
    animalFractionImproved;


AllCSMAnimalTable = table( ...
    Animal, ...
    T0Reaches, ...
    T1Reaches, ...
    T2Reaches, ...
    T3Reaches, ...
    TotalReaches, ...
    IncludedInStatistics, ...
    PredictedMSE, ...
    ActualMSE, ...
    PredictedRMSE, ...
    ActualRMSE, ...
    MeanImprovement, ...
    FractionImproved);


disp(AllCSMAnimalTable);


%% ========================================================================
% STORE COMBINED RESULTS
% ========================================================================

AllCSMEndpointResults = struct;


%% Settings

AllCSMEndpointResults.settings.subtypes = ...
    combinedSubtypes;

AllCSMEndpointResults.settings.positionScale = ...
    positionScale;

AllCSMEndpointResults.settings.positionUnit = ...
    positionUnit;

AllCSMEndpointResults.settings.minimumReachesPerAnimal = ...
    minimumCombinedReachesPerAnimal;

AllCSMEndpointResults.settings.targetPosition = ...
    targetPosition;


%% Reach-level results

AllCSMEndpointResults.reach.predictedEndpoint = ...
    allPredictedEndpoints;

AllCSMEndpointResults.reach.actualEndpoint = ...
    allActualEndpoints;

AllCSMEndpointResults.reach.predictedSquaredError = ...
    predictedSquaredError;

AllCSMEndpointResults.reach.actualSquaredError = ...
    actualSquaredError;

AllCSMEndpointResults.reach.predictedEuclideanError = ...
    predictedEuclideanError;

AllCSMEndpointResults.reach.actualEuclideanError = ...
    actualEuclideanError;

AllCSMEndpointResults.reach.improvement = ...
    endpointErrorImprovement;

AllCSMEndpointResults.reach.squaredErrorImprovement = ...
    squaredErrorImprovement;

AllCSMEndpointResults.reach.animalID = ...
    allAnimalIDs;

AllCSMEndpointResults.reach.subtypeID = ...
    allSubtypeIDs;

AllCSMEndpointResults.reach.subtypeLabel = ...
    allSubtypeLabels;

AllCSMEndpointResults.reach.originalReachIndex = ...
    allOriginalReachIndices;


%% Animal-level results

AllCSMEndpointResults.animal.ID = ...
    uniqueAnimalIDs;

AllCSMEndpointResults.animal.reachCount = ...
    animalCombinedReachCount;

AllCSMEndpointResults.animal.subtypeReachCounts = ...
    animalSubtypeReachCounts;

AllCSMEndpointResults.animal.predictedMSE = ...
    animalPredictedMSE;

AllCSMEndpointResults.animal.actualMSE = ...
    animalActualMSE;

AllCSMEndpointResults.animal.predictedRMSE = ...
    animalPredictedRMSE;

AllCSMEndpointResults.animal.actualRMSE = ...
    animalActualRMSE;

AllCSMEndpointResults.animal.meanPredictedError = ...
    animalMeanPredictedError;

AllCSMEndpointResults.animal.meanActualError = ...
    animalMeanActualError;

AllCSMEndpointResults.animal.meanImprovement = ...
    animalMeanImprovement;

AllCSMEndpointResults.animal.fractionImproved = ...
    animalFractionImproved;

AllCSMEndpointResults.animal.validForMSE = ...
    validAnimalsMSE;


%% Group results

AllCSMEndpointResults.group.predictedMSE = ...
    groupPredictedMSE;

AllCSMEndpointResults.group.actualMSE = ...
    groupActualMSE;

AllCSMEndpointResults.group.predictedRMSE = ...
    groupPredictedRMSE;

AllCSMEndpointResults.group.actualRMSE = ...
    groupActualRMSE;

AllCSMEndpointResults.group.predictedEuclideanError = ...
    groupPredictedEuclideanError;

AllCSMEndpointResults.group.actualEuclideanError = ...
    groupActualEuclideanError;

AllCSMEndpointResults.group.meanImprovement = ...
    groupEndpointImprovement;

AllCSMEndpointResults.group.fractionImproved = ...
    groupFractionImproved;


%% Statistics

AllCSMEndpointResults.statistics.pMSE = ...
    pCombinedMSE;

AllCSMEndpointResults.statistics.pEuclideanError = ...
    pCombinedEuclideanError;

AllCSMEndpointResults.statistics.pImprovement = ...
    pCombinedImprovement;

AllCSMEndpointResults.statistics.pFractionImproved = ...
    pCombinedFractionImproved;


%% Tables

AllCSMEndpointResults.animalTable = ...
    AllCSMAnimalTable;


fprintf('\nCombined analysis complete.\n');
fprintf('Results stored in AllCSMEndpointResults.\n');


%% ========================================================================
% CREATE SUMMARY TABLE
% ========================================================================

Subtype = subtypes(:);

NAnimals = numberAnimalsBySubtype;
NReaches = numberReachesBySubtype;

PredictedMSE = predictedMSEBySubtype;
ActualMSE = actualMSEBySubtype;

PredictedRMSE = predictedRMSEBySubtype;
ActualRMSE = actualRMSEBySubtype;

pMSE = pMSEBySubtype;
qMSE = qMSEBySubtype;


EndpointTargetErrorTable = table( ...
    Subtype, ...
    NAnimals, ...
    NReaches, ...
    PredictedMSE, ...
    ActualMSE, ...
    PredictedRMSE, ...
    ActualRMSE, ...
    pMSE, ...
    qMSE);


disp(EndpointTargetErrorTable);


%% Store with the primary results structure

DecelerationEndpointResults.EndpointTargetStats = ...
    SubtypeEndpointStats;

DecelerationEndpointResults.EndpointTargetErrorTable = ...
    EndpointTargetErrorTable;

%% ========================================================================
% LOCAL FUNCTION: BENJAMINI-HOCHBERG FDR
% ========================================================================

function adjustedP = benjaminiHochbergLocal(pValues)


    originalSize = size(pValues);

    pValues = pValues(:);

    adjustedP = nan(size(pValues));


    valid = isfinite(pValues);

    validP = pValues(valid);


    if isempty(validP)

        adjustedP = reshape(adjustedP,originalSize);

        return

    end


    [sortedP,sortOrder] = sort(validP);

    numberTests = numel(sortedP);


    sortedAdjusted = ...
        sortedP .* numberTests ./ (1:numberTests)';


    for idx = numberTests-1:-1:1

        sortedAdjusted(idx) = min( ...
            sortedAdjusted(idx), ...
            sortedAdjusted(idx+1));

    end


    sortedAdjusted = min(sortedAdjusted,1);


    unsortedAdjusted = nan(numberTests,1);

    unsortedAdjusted(sortOrder) = ...
        sortedAdjusted;


    adjustedP(valid) = ...
        unsortedAdjusted;


    adjustedP = ...
        reshape(adjustedP,originalSize);

end

%% 
%% ========================================================================
% PREDICTED VERSUS ACTUAL OUTWARD ENDPOINT POSITION
% ALL CSM SUBTYPES COMBINED
%
% Predicted:
%   Maximum-outward endpoint predicted from the pre-CSM velocity and
%   deceleration model.
%
% Actual:
%   Measured paw position at maximum outward extension.
%
% For each animal:
%
%   1. Pool all valid T0-T3 reaches.
%   2. Extract the outward coordinate of the predicted endpoint.
%   3. Extract the outward coordinate of the actual endpoint.
%   4. Average predicted and actual values within the animal.
%   5. Plot paired animal means.
%   6. Plot the grand mean +/- SEM.
%   7. Perform a paired t-test.
%   8. Draw a horizontal line at the target outward position.
%
% REQUIRED:
%
%   DecelerationEndpointResults.T0
%   DecelerationEndpointResults.T1
%   DecelerationEndpointResults.T2
%   DecelerationEndpointResults.T3
% ========================================================================


%% ------------------------------------------------------------------------
% SETTINGS
% -------------------------------------------------------------------------

combinedSubtypes = {'T0','T1','T2','T3'};

% Endpoint arrays contain the dimensions selected in analysisDimensions.
%
% With:
%
%   analysisDimensions = [1 2]
%
% endpoint column 1 = x
% endpoint column 2 = y
%
% Set to 1 to compare outward x position.

outwardEndpointColumn = 1;


%% Units

% Original data are in decimeters.
%
% Centimeters:
positionScale = 10;
positionUnit = 'cm';

% Decimeters:
% positionScale = 1;
% positionUnit = 'dm';


%% Animal inclusion

minimumCombinedReachesPerAnimal = 5;


%% Target outward coordinate

if ~exist('targetPosition','var') || ...
        ~exist('analysisDimensions','var')

    error(['targetPosition and analysisDimensions must already exist ' ...
        'in the workspace.']);

end

target2D = ...
    targetPosition(analysisDimensions);

targetOutwardPosition = ...
    target2D(outwardEndpointColumn) .* positionScale;


%% ========================================================================
% COMBINE T0-T3 REACHES
% ========================================================================

allPredictedOutward = [];
allActualOutward = [];

allAnimalIDs = [];
allSubtypeIDs = [];
allOriginalReachIndices = [];


for subtypeIdx = 1:numel(combinedSubtypes)

    subtype = combinedSubtypes{subtypeIdx};


    if ~isfield(DecelerationEndpointResults,subtype)

        warning('%s is missing from DecelerationEndpointResults.', ...
            subtype);

        continue

    end


    R = DecelerationEndpointResults.(subtype);


    if ~isfield(R,'reach') || ...
            ~isfield(R.reach,'predictedEndpoint') || ...
            ~isfield(R.reach,'actualEndpoint') || ...
            ~isfield(R.reach,'animalID')

        warning('%s does not contain the required reach-level fields.', ...
            subtype);

        continue

    end


    predictedEndpoints = ...
        R.reach.predictedEndpoint;

    actualEndpoints = ...
        R.reach.actualEndpoint;

    animalIDs = ...
        R.reach.animalID;


    if isfield(R.reach,'originalReachIndex')

        originalReachIndices = ...
            R.reach.originalReachIndex;

    else

        originalReachIndices = ...
            (1:size(predictedEndpoints,1))';

    end


    %% Retain valid paired endpoint values

    validPairs = ...
        all(isfinite(predictedEndpoints),2) & ...
        all(isfinite(actualEndpoints),2) & ...
        isfinite(animalIDs);


    predictedEndpoints = ...
        predictedEndpoints(validPairs,:);

    actualEndpoints = ...
        actualEndpoints(validPairs,:);

    animalIDs = ...
        animalIDs(validPairs);

    originalReachIndices = ...
        originalReachIndices(validPairs);


    predictedOutward = ...
        predictedEndpoints(:,outwardEndpointColumn);

    actualOutward = ...
        actualEndpoints(:,outwardEndpointColumn);


    numberSubtypeReaches = ...
        numel(predictedOutward);


    allPredictedOutward = [ ...
        allPredictedOutward; ...
        predictedOutward];

    allActualOutward = [ ...
        allActualOutward; ...
        actualOutward];

    allAnimalIDs = [ ...
        allAnimalIDs; ...
        animalIDs];

    allSubtypeIDs = [ ...
        allSubtypeIDs; ...
        repmat(subtypeIdx,numberSubtypeReaches,1)];

    allOriginalReachIndices = [ ...
        allOriginalReachIndices; ...
        originalReachIndices];


    fprintf('%s: added %d valid reaches\n', ...
        subtype,numberSubtypeReaches);

end


if isempty(allPredictedOutward)

    error('No valid predicted and actual endpoints were found.');

end


fprintf('\nTotal combined reaches: %d\n', ...
    numel(allPredictedOutward));


%% ========================================================================
% ANIMAL-LEVEL MEANS
% ========================================================================

uniqueAnimalIDs = ...
    unique(allAnimalIDs);

numberCandidateAnimals = ...
    numel(uniqueAnimalIDs);


animalPredictedMean = ...
    nan(numberCandidateAnimals,1);

animalActualMean = ...
    nan(numberCandidateAnimals,1);

animalReachCount = ...
    zeros(numberCandidateAnimals,1);

animalSubtypeReachCount = ...
    zeros(numberCandidateAnimals,numel(combinedSubtypes));


for animalIdx = 1:numberCandidateAnimals

    thisAnimalID = ...
        uniqueAnimalIDs(animalIdx);

    animalMask = ...
        allAnimalIDs == thisAnimalID;


    animalReachCount(animalIdx) = ...
        sum(animalMask);


    for subtypeIdx = 1:numel(combinedSubtypes)

        animalSubtypeReachCount(animalIdx,subtypeIdx) = ...
            sum( ...
                animalMask & ...
                allSubtypeIDs == subtypeIdx);

    end


    if animalReachCount(animalIdx) < ...
            minimumCombinedReachesPerAnimal

        continue

    end


    animalPredictedMean(animalIdx) = ...
        mean( ...
            allPredictedOutward(animalMask), ...
            'omitnan');


    animalActualMean(animalIdx) = ...
        mean( ...
            allActualOutward(animalMask), ...
            'omitnan');

end


%% Retain paired animals

validAnimals = ...
    isfinite(animalPredictedMean) & ...
    isfinite(animalActualMean);


pairedPredicted = ...
    animalPredictedMean(validAnimals) .* positionScale;

pairedActual = ...
    animalActualMean(validAnimals) .* positionScale;

pairedAnimalIDs = ...
    uniqueAnimalIDs(validAnimals);

pairedReachCount = ...
    animalReachCount(validAnimals);

pairedSubtypeReachCount = ...
    animalSubtypeReachCount(validAnimals,:);


numberPairedAnimals = ...
    sum(validAnimals);


if numberPairedAnimals < 2

    error(['Fewer than two animals had sufficient paired predicted ' ...
        'and actual endpoint data.']);

end


%% ========================================================================
% PAIRED T-TEST
%
% The default paired t-test is two-sided:
%
%   H0: mean predicted outward position =
%       mean actual outward position
% ========================================================================

[hPaired,pPaired,confidenceInterval,tTestStats] = ...
    ttest(pairedPredicted,pairedActual);


%% Difference defined as actual - predicted

animalDifference = ...
    pairedActual - pairedPredicted;


meanDifference = ...
    mean(animalDifference,'omitnan');

semDifference = ...
    std(animalDifference,0,'omitnan') ./ ...
    sqrt(numberPairedAnimals);


%% Group means and SEM

grandMeanPredicted = ...
    mean(pairedPredicted,'omitnan');

grandMeanActual = ...
    mean(pairedActual,'omitnan');


semPredicted = ...
    std(pairedPredicted,0,'omitnan') ./ ...
    sqrt(numberPairedAnimals);

semActual = ...
    std(pairedActual,0,'omitnan') ./ ...
    sqrt(numberPairedAnimals);


%% ========================================================================
% PRINT RESULTS
% ========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('PREDICTED VERSUS ACTUAL OUTWARD ENDPOINT\n');
fprintf('ALL CSM SUBTYPES COMBINED\n');
fprintf('============================================================\n');

fprintf('Paired animals: %d\n',numberPairedAnimals);

fprintf('Total reaches from included animals: %d\n', ...
    sum(pairedReachCount));


fprintf('\nPredicted outward endpoint: %.5f +/- %.5f %s\n', ...
    grandMeanPredicted, ...
    semPredicted, ...
    positionUnit);


fprintf('Actual outward endpoint: %.5f +/- %.5f %s\n', ...
    grandMeanActual, ...
    semActual, ...
    positionUnit);


fprintf('Actual - predicted difference: %.5f +/- %.5f %s\n', ...
    meanDifference, ...
    semDifference, ...
    positionUnit);


fprintf('Paired t(%d) = %.4f, p = %.6g\n', ...
    tTestStats.df, ...
    tTestStats.tstat, ...
    pPaired);


fprintf(['95%% CI of predicted - actual difference: ' ...
    '[%.5f, %.5f] %s\n'], ...
    confidenceInterval(1), ...
    confidenceInterval(2), ...
    positionUnit);


fprintf('Target outward position: %.5f %s\n', ...
    targetOutwardPosition, ...
    positionUnit);


%% ========================================================================
% PLOT
% ========================================================================

figurePredictedVsActualOutward = figure( ...
    'Name','Predicted versus actual outward endpoint', ...
    'Color','w', ...
    'Position',[150 100 650 650]);


ax = axes(figurePredictedVsActualOutward);

hold(ax,'on');


%% ------------------------------------------------------------------------
% Paired animal means
% -------------------------------------------------------------------------

for animalIdx = 1:numberPairedAnimals

    plot( ...
        ax, ...
        [1 2], ...
        [ ...
            pairedPredicted(animalIdx), ...
            pairedActual(animalIdx)], ...
        '-o', ...
        'Color',[0.65 0.65 0.65], ...
        'MarkerFaceColor','w', ...
        'MarkerEdgeColor',[0.35 0.35 0.35], ...
        'LineWidth',1.3, ...
        'MarkerSize',7);

end


%% ------------------------------------------------------------------------
% Grand means +/- SEM
% -------------------------------------------------------------------------

grandMeanHandle = errorbar( ...
    ax, ...
    [1 2], ...
    [ ...
        grandMeanPredicted, ...
        grandMeanActual], ...
    [ ...
        semPredicted, ...
        semActual], ...
    'ko-', ...
    'MarkerFaceColor','k', ...
    'MarkerSize',10, ...
    'LineWidth',2.6, ...
    'CapSize',12);


%% ------------------------------------------------------------------------
% Target line
% -------------------------------------------------------------------------

targetLineHandle = yline( ...
    ax, ...
    targetOutwardPosition, ...
    'k--', ...
    'Target', ...
    'LineWidth',1.8, ...
    'LabelHorizontalAlignment','right', ...
    'LabelVerticalAlignment','bottom');


%% ------------------------------------------------------------------------
% Formatting
% -------------------------------------------------------------------------

xlim(ax,[0.65 2.35]);

xticks(ax,[1 2]);

xticklabels(ax,{ ...
    'Predicted', ...
    'Actual'});


ylabel(ax,sprintf( ...
    'Maximum outward position (%s)', ...
    positionUnit));


title(ax,{ ...
    sprintf('All CSM subtypes combined, N = %d animals', ...
        numberPairedAnimals), ...
    sprintf('Paired t(%d) = %.2f, p = %.3g', ...
        tTestStats.df, ...
        tTestStats.tstat, ...
        pPaired)});


box(ax,'off');


set( ...
    ax, ...
    'FontSize',12, ...
    'LineWidth',1.3, ...
    'TickDir','out');


%% Set y-limits so all data and the target are visible

allPlotValues = [ ...
    pairedPredicted; ...
    pairedActual; ...
    targetOutwardPosition];


plotMinimum = ...
    min(allPlotValues);

plotMaximum = ...
    max(allPlotValues);

plotRange = ...
    plotMaximum - plotMinimum;


if plotRange == 0

    plotRange = max(abs(allPlotValues))*0.1;

end

if plotRange == 0

    plotRange = 1;

end


ylim(ax,[ ...
    plotMinimum - 0.12*plotRange, ...
    plotMaximum + 0.12*plotRange]);


%% ========================================================================
% ANIMAL TABLE
% ========================================================================

if exist('animalNames','var') && ...
        numel(animalNames) >= max(pairedAnimalIDs)

    Animal = ...
        animalNames(pairedAnimalIDs);

else

    Animal = ...
        compose('Animal_%d',pairedAnimalIDs);

end


T0ReachCount = ...
    pairedSubtypeReachCount(:,1);

T1ReachCount = ...
    pairedSubtypeReachCount(:,2);

T2ReachCount = ...
    pairedSubtypeReachCount(:,3);

T3ReachCount = ...
    pairedSubtypeReachCount(:,4);

TotalReachCount = ...
    pairedReachCount;

PredictedMeanOutwardPosition = ...
    pairedPredicted;

ActualMeanOutwardPosition = ...
    pairedActual;

ActualMinusPredicted = ...
    animalDifference;


PredictedVsActualOutwardTable = table( ...
    Animal, ...
    T0ReachCount, ...
    T1ReachCount, ...
    T2ReachCount, ...
    T3ReachCount, ...
    TotalReachCount, ...
    PredictedMeanOutwardPosition, ...
    ActualMeanOutwardPosition, ...
    ActualMinusPredicted);


disp(PredictedVsActualOutwardTable);


%% ========================================================================
% STORE RESULTS
% ========================================================================

PredictedVsActualOutwardResults = struct;


PredictedVsActualOutwardResults.settings.subtypes = ...
    combinedSubtypes;

PredictedVsActualOutwardResults.settings.outwardEndpointColumn = ...
    outwardEndpointColumn;

PredictedVsActualOutwardResults.settings.positionScale = ...
    positionScale;

PredictedVsActualOutwardResults.settings.positionUnit = ...
    positionUnit;

PredictedVsActualOutwardResults.settings.targetOutwardPosition = ...
    targetOutwardPosition;

PredictedVsActualOutwardResults.settings.minimumReachesPerAnimal = ...
    minimumCombinedReachesPerAnimal;


PredictedVsActualOutwardResults.animal.ID = ...
    pairedAnimalIDs;

PredictedVsActualOutwardResults.animal.predictedMean = ...
    pairedPredicted;

PredictedVsActualOutwardResults.animal.actualMean = ...
    pairedActual;

PredictedVsActualOutwardResults.animal.actualMinusPredicted = ...
    animalDifference;

PredictedVsActualOutwardResults.animal.reachCount = ...
    pairedReachCount;

PredictedVsActualOutwardResults.animal.subtypeReachCount = ...
    pairedSubtypeReachCount;


PredictedVsActualOutwardResults.group.predictedMean = ...
    grandMeanPredicted;

PredictedVsActualOutwardResults.group.actualMean = ...
    grandMeanActual;

PredictedVsActualOutwardResults.group.predictedSEM = ...
    semPredicted;

PredictedVsActualOutwardResults.group.actualSEM = ...
    semActual;

PredictedVsActualOutwardResults.group.meanDifference = ...
    meanDifference;

PredictedVsActualOutwardResults.group.semDifference = ...
    semDifference;


PredictedVsActualOutwardResults.statistics.h = ...
    hPaired;

PredictedVsActualOutwardResults.statistics.p = ...
    pPaired;

PredictedVsActualOutwardResults.statistics.t = ...
    tTestStats.tstat;

PredictedVsActualOutwardResults.statistics.df = ...
    tTestStats.df;

PredictedVsActualOutwardResults.statistics.confidenceInterval = ...
    confidenceInterval;


PredictedVsActualOutwardResults.animalTable = ...
    PredictedVsActualOutwardTable;


fprintf('\nAnalysis complete.\n');
fprintf('Results stored in PredictedVsActualOutwardResults.\n');

%% ========================================================================
% PREDICTED VERSUS ACTUAL SIGNED DIFFERENCE FROM TARGET
% ALL CSM SUBTYPES COMBINED
%
% Signed target difference:
%
%   endpoint position - target position
%
% Interpretation:
%
%   0        = endpoint is at the target
%   positive = endpoint extends beyond the target
%   negative = endpoint falls short of the target
%
% Uses variables generated by the previous predicted-versus-actual analysis:
%
%   pairedPredicted
%   pairedActual
%   targetOutwardPosition
%   numberPairedAnimals
%   pairedAnimalIDs
%   pairedReachCount
%   positionUnit
% ========================================================================


%% Calculate signed difference from target

predictedDifferenceFromTarget = ...
    pairedPredicted - targetOutwardPosition;

actualDifferenceFromTarget = ...
    pairedActual - targetOutwardPosition;


%% ========================================================================
% PAIRED T-TEST
%
% Two-sided null hypothesis:
%
%   predicted signed target difference =
%   actual signed target difference
% ========================================================================

[hDifference,pDifference,ciDifference,differenceTTestStats] = ...
    ttest( ...
        predictedDifferenceFromTarget, ...
        actualDifferenceFromTarget);


%% Difference between conditions
%
% Positive values mean the actual endpoint is farther outward than predicted.
% Negative values mean the actual endpoint is less outward than predicted.

animalActualMinusPredicted = ...
    actualDifferenceFromTarget - ...
    predictedDifferenceFromTarget;


meanActualMinusPredicted = ...
    mean(animalActualMinusPredicted,'omitnan');

semActualMinusPredicted = ...
    std(animalActualMinusPredicted,0,'omitnan') ./ ...
    sqrt(numberPairedAnimals);


%% Grand means and SEM

grandMeanPredictedDifference = ...
    mean(predictedDifferenceFromTarget,'omitnan');

grandMeanActualDifference = ...
    mean(actualDifferenceFromTarget,'omitnan');


semPredictedDifference = ...
    std(predictedDifferenceFromTarget,0,'omitnan') ./ ...
    sqrt(numberPairedAnimals);

semActualDifference = ...
    std(actualDifferenceFromTarget,0,'omitnan') ./ ...
    sqrt(numberPairedAnimals);


%% ========================================================================
% PRINT RESULTS
% ========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('PREDICTED VERSUS ACTUAL SIGNED DIFFERENCE FROM TARGET\n');
fprintf('============================================================\n');

fprintf('Animals: %d\n',numberPairedAnimals);


fprintf('\nPredicted difference from target: %.5f +/- %.5f %s\n', ...
    grandMeanPredictedDifference, ...
    semPredictedDifference, ...
    positionUnit);


fprintf('Actual difference from target: %.5f +/- %.5f %s\n', ...
    grandMeanActualDifference, ...
    semActualDifference, ...
    positionUnit);


fprintf('Actual - predicted difference: %.5f +/- %.5f %s\n', ...
    meanActualMinusPredicted, ...
    semActualMinusPredicted, ...
    positionUnit);


fprintf('Paired t(%d) = %.4f, p = %.6g\n', ...
    differenceTTestStats.df, ...
    differenceTTestStats.tstat, ...
    pDifference);


fprintf(['95%% CI of predicted - actual difference: ' ...
    '[%.5f, %.5f] %s\n'], ...
    ciDifference(1), ...
    ciDifference(2), ...
    positionUnit);


%% ========================================================================
% NARROW PAIRED PLOT
% ========================================================================

% figureSignedDifference = figure( ...
%     'Name','Predicted versus actual signed target difference', ...
%     'Color','w', ...
%     'Position',[180 100 360 650]);


figureSignedDifference = figure( ...
    'Name','Predicted versus actual signed target difference', ...
    'Color','w', ...
    'Position',[200 200 420 900]);



axDifference = axes(figureSignedDifference);

hold(axDifference,'on');


%% Paired animal values

for animalIdx = 1:numberPairedAnimals

    plot( ...
        axDifference, ...
        [1 2], ...
        [ ...
            predictedDifferenceFromTarget(animalIdx), ...
            actualDifferenceFromTarget(animalIdx)], ...
        '-o', ...
        'Color',[0.65 0.65 0.65], ...
        'MarkerFaceColor','w', ...
        'MarkerEdgeColor',[0.35 0.35 0.35], ...
        'LineWidth',1.3, ...
        'MarkerSize',7);

end


%% Grand mean +/- SEM

errorbar( ...
    axDifference, ...
    [1 2], ...
    [ ...
        grandMeanPredictedDifference, ...
        grandMeanActualDifference], ...
    [ ...
        semPredictedDifference, ...
        semActualDifference], ...
    'ko-', ...
    'MarkerFaceColor','k', ...
    'MarkerSize',10, ...
    'LineWidth',2.6, ...
    'CapSize',10);


%% Zero line marks the target

yline( ...
    axDifference, ...
    0, ...
    'k--', ...
    'Target', ...
    'LineWidth',1.6, ...
    'LabelHorizontalAlignment','right', ...
    'LabelVerticalAlignment','bottom');


%% Formatting

xlim(axDifference,[0.82 2.18]);

xticks(axDifference,[1 2]);

xticklabels(axDifference,{ ...
    'Predicted', ...
    'Actual'});


ylabel(axDifference,sprintf( ...
    'Signed outward difference from target (%s)', ...
    positionUnit));


title(axDifference,{ ...
    sprintf('All CSMs, N = %d animals',numberPairedAnimals), ...
    sprintf('Paired t(%d) = %.2f, p = %.3g', ...
        differenceTTestStats.df, ...
        differenceTTestStats.tstat, ...
        pDifference)});


box(axDifference,'off');


set( ...
    axDifference, ...
    'FontSize',12, ...
    'LineWidth',1.3, ...
    'TickDir','out');


%% Symmetric y-axis around zero

maximumAbsoluteDifference = max(abs([ ...
    predictedDifferenceFromTarget; ...
    actualDifferenceFromTarget]));


if ~isfinite(maximumAbsoluteDifference) || ...
        maximumAbsoluteDifference == 0

    maximumAbsoluteDifference = 1;

end


ylim(axDifference,[ ...
    -1.15 * maximumAbsoluteDifference, ...
     1.15 * maximumAbsoluteDifference]);


%% ========================================================================
% TABLE
% ========================================================================

if exist('animalNames','var') && ...
        numel(animalNames) >= max(pairedAnimalIDs)

    Animal = animalNames(pairedAnimalIDs);

else

    Animal = compose('Animal_%d',pairedAnimalIDs);

end


PredictedDifferenceFromTarget = ...
    predictedDifferenceFromTarget;

ActualDifferenceFromTarget = ...
    actualDifferenceFromTarget;

ActualMinusPredicted = ...
    animalActualMinusPredicted;

TotalReachCount = ...
    pairedReachCount;


SignedDifferenceFromTargetTable = table( ...
    Animal, ...
    TotalReachCount, ...
    PredictedDifferenceFromTarget, ...
    ActualDifferenceFromTarget, ...
    ActualMinusPredicted);


disp(SignedDifferenceFromTargetTable);


%% ========================================================================
% STORE RESULTS
% ========================================================================

SignedDifferenceFromTargetResults = struct;


SignedDifferenceFromTargetResults.settings.positionUnit = ...
    positionUnit;

SignedDifferenceFromTargetResults.settings.targetOutwardPosition = ...
    targetOutwardPosition;


SignedDifferenceFromTargetResults.animal.ID = ...
    pairedAnimalIDs;

SignedDifferenceFromTargetResults.animal.predictedDifference = ...
    predictedDifferenceFromTarget;

SignedDifferenceFromTargetResults.animal.actualDifference = ...
    actualDifferenceFromTarget;

SignedDifferenceFromTargetResults.animal.actualMinusPredicted = ...
    animalActualMinusPredicted;

SignedDifferenceFromTargetResults.animal.reachCount = ...
    pairedReachCount;


SignedDifferenceFromTargetResults.group.predictedMean = ...
    grandMeanPredictedDifference;

SignedDifferenceFromTargetResults.group.actualMean = ...
    grandMeanActualDifference;

SignedDifferenceFromTargetResults.group.predictedSEM = ...
    semPredictedDifference;

SignedDifferenceFromTargetResults.group.actualSEM = ...
    semActualDifference;

SignedDifferenceFromTargetResults.group.meanActualMinusPredicted = ...
    meanActualMinusPredicted;

SignedDifferenceFromTargetResults.group.semActualMinusPredicted = ...
    semActualMinusPredicted;


SignedDifferenceFromTargetResults.statistics.h = ...
    hDifference;

SignedDifferenceFromTargetResults.statistics.p = ...
    pDifference;

SignedDifferenceFromTargetResults.statistics.t = ...
    differenceTTestStats.tstat;

SignedDifferenceFromTargetResults.statistics.df = ...
    differenceTTestStats.df;

SignedDifferenceFromTargetResults.statistics.confidenceInterval = ...
    ciDifference;


SignedDifferenceFromTargetResults.animalTable = ...
    SignedDifferenceFromTargetTable;


fprintf('\nAnalysis complete.\n');
fprintf('Results stored in SignedDifferenceFromTargetResults.\n');

%% other option 

%% ========================================================================
% REACH-LEVEL PREDICTED VERSUS ACTUAL SIGNED DIFFERENCE FROM TARGET
% ALL CSM SUBTYPES COMBINED
%
% This code:
%
%   1. Combines valid reaches from T0, T1, T2, and T3.
%   2. Extracts predicted and actual outward endpoint position.
%   3. Calculates signed difference from the target:
%
%        endpoint outward position - target outward position
%
%   4. Plots every reach as a connected predicted-actual pair.
%   5. Plots predicted values in dark gray.
%   6. Plots actual values in purple.
%   7. Draws paired reach lines in light gray behind the points.
%   8. Plots grand mean +/- SEM over all reaches.
%   9. Performs a reach-level paired t-test.
%
% Interpretation:
%
%   0        = endpoint is at target
%   negative = endpoint falls short of target
%   positive = endpoint extends beyond target
%
% REQUIRED WORKSPACE VARIABLES:
%
%   DecelerationEndpointResults
%   targetPosition
%   analysisDimensions
%
% Expected fields:
%
%   DecelerationEndpointResults.T0.reach.predictedEndpoint
%   DecelerationEndpointResults.T0.reach.actualEndpoint
%   DecelerationEndpointResults.T0.reach.animalID
%
% and equivalent fields for T1-T3.
% ========================================================================


%% ========================================================================
% SETTINGS
% ========================================================================

combinedSubtypes = {'T0','T1','T2','T3'};


%% Endpoint coordinate to analyze
%
% Endpoint arrays contain the coordinates selected by analysisDimensions.
%
% With:
%
%   analysisDimensions = [1 2]
%
% endpoint column 1 = x
% endpoint column 2 = y
%
% Set this to 1 for outward x position.

outwardEndpointColumn = 1;


%% Units
%
% Original positions are in decimeters.
%
% For centimeters:

positionScale = 10;
positionUnit = 'cm';

% For decimeters instead, use:
%
% positionScale = 1;
% positionUnit = 'dm';


%% Plot appearance

predictedPointColor = [0.25 0.25 0.25];

actualPointColor = [0.5781 0 0.8242];

pairedLineColor = [0.82 0.82 0.82];

predictedPointSize = 18;

actualPointSize = 18;

pairedLineWidth = 0.40;

individualPointAlpha = 0.32;

grandMeanMarkerSize = 11;

grandMeanLineWidth = 2.8;


%% Figure size
%
% Increase the third value to make the figure wider.

%figurePosition = [180 100 420 900];

figurePosition = [180 100 400 300];


%% Optional manually specified y-limits
%
% Leave as [NaN NaN] for automatic symmetric limits around zero.
%
% Example:
%
% manualYLimits = [-1.5 1.5];

manualYLimits = [NaN NaN];


%% ========================================================================
% VALIDATE REQUIRED VARIABLES
% ========================================================================

if ~exist('DecelerationEndpointResults','var')

    error('DecelerationEndpointResults is missing from the workspace.');

end


if ~exist('targetPosition','var') || ...
        numel(targetPosition) < 2

    error('targetPosition is missing or invalid.');

end


if ~exist('analysisDimensions','var') || ...
        isempty(analysisDimensions)

    error('analysisDimensions is missing from the workspace.');

end


if outwardEndpointColumn > numel(analysisDimensions)

    error(['outwardEndpointColumn exceeds the number of coordinates ' ...
        'stored in the endpoint arrays.']);

end


%% ========================================================================
% DETERMINE TARGET OUTWARD POSITION
% ========================================================================

target2D = ...
    targetPosition(analysisDimensions);


targetOutwardPositionOriginalUnits = ...
    target2D(outwardEndpointColumn);


targetOutwardPositionForPlot = ...
    targetOutwardPositionOriginalUnits .* positionScale;


fprintf('\nTarget outward position: %.5f %s\n', ...
    targetOutwardPositionForPlot, ...
    positionUnit);


%% ========================================================================
% COMBINE VALID REACHES FROM T0-T3
% ========================================================================

allPredictedOutward = [];

allActualOutward = [];

allAnimalIDs = [];

allSubtypeIDs = [];

allSubtypeLabels = {};

allOriginalReachIndices = [];


for subtypeIdx = 1:numel(combinedSubtypes)

    subtype = combinedSubtypes{subtypeIdx};


    %% Confirm subtype exists

    if ~isfield(DecelerationEndpointResults,subtype)

        warning('%s is missing from DecelerationEndpointResults.', ...
            subtype);

        continue

    end


    R = DecelerationEndpointResults.(subtype);


    %% Confirm required reach fields exist

    if ~isfield(R,'reach') || ...
            ~isfield(R.reach,'predictedEndpoint') || ...
            ~isfield(R.reach,'actualEndpoint')

        warning('%s does not contain predicted and actual endpoints.', ...
            subtype);

        continue

    end


    predictedEndpoints = ...
        R.reach.predictedEndpoint;

    actualEndpoints = ...
        R.reach.actualEndpoint;


    %% Animal identity

    if isfield(R.reach,'animalID')

        animalIDs = ...
            R.reach.animalID;

    else

        animalIDs = ...
            nan(size(predictedEndpoints,1),1);

    end


    %% Original reach index

    if isfield(R.reach,'originalReachIndex')

        originalReachIndices = ...
            R.reach.originalReachIndex;

    else

        originalReachIndices = ...
            (1:size(predictedEndpoints,1))';

    end


    %% Ensure arrays have matching lengths

    numberRows = min([ ...
        size(predictedEndpoints,1), ...
        size(actualEndpoints,1), ...
        numel(animalIDs), ...
        numel(originalReachIndices)]);


    predictedEndpoints = ...
        predictedEndpoints(1:numberRows,:);

    actualEndpoints = ...
        actualEndpoints(1:numberRows,:);

    animalIDs = ...
        animalIDs(1:numberRows);

    originalReachIndices = ...
        originalReachIndices(1:numberRows);


    %% Retain valid matched endpoint pairs

    validPairs = ...
        all(isfinite(predictedEndpoints),2) & ...
        all(isfinite(actualEndpoints),2);


    predictedEndpoints = ...
        predictedEndpoints(validPairs,:);

    actualEndpoints = ...
        actualEndpoints(validPairs,:);

    animalIDs = ...
        animalIDs(validPairs);

    originalReachIndices = ...
        originalReachIndices(validPairs);


    %% Extract outward coordinate

    predictedOutward = ...
        predictedEndpoints(:,outwardEndpointColumn);

    actualOutward = ...
        actualEndpoints(:,outwardEndpointColumn);


    validOutwardPairs = ...
        isfinite(predictedOutward) & ...
        isfinite(actualOutward);


    predictedOutward = ...
        predictedOutward(validOutwardPairs);

    actualOutward = ...
        actualOutward(validOutwardPairs);

    animalIDs = ...
        animalIDs(validOutwardPairs);

    originalReachIndices = ...
        originalReachIndices(validOutwardPairs);


    numberSubtypeReaches = ...
        numel(predictedOutward);


    %% Add to combined arrays

    allPredictedOutward = [ ...
        allPredictedOutward; ...
        predictedOutward];


    allActualOutward = [ ...
        allActualOutward; ...
        actualOutward];


    allAnimalIDs = [ ...
        allAnimalIDs; ...
        animalIDs];


    allSubtypeIDs = [ ...
        allSubtypeIDs; ...
        repmat(subtypeIdx,numberSubtypeReaches,1)];


    allSubtypeLabels = [ ...
        allSubtypeLabels; ...
        repmat({subtype},numberSubtypeReaches,1)];


    allOriginalReachIndices = [ ...
        allOriginalReachIndices; ...
        originalReachIndices];


    fprintf('%s: added %d valid reaches\n', ...
        subtype, ...
        numberSubtypeReaches);

end


%% Confirm data exist

numberValidReaches = ...
    numel(allPredictedOutward);


if numberValidReaches < 2

    error('Fewer than two valid paired reaches were found.');

end


fprintf('\nTotal valid paired reaches: %d\n', ...
    numberValidReaches);


%% ========================================================================
% CONVERT TO DISPLAY UNITS
% ========================================================================

predictedOutwardForPlot = ...
    allPredictedOutward .* positionScale;


actualOutwardForPlot = ...
    allActualOutward .* positionScale;


%% ========================================================================
% CALCULATE SIGNED DIFFERENCE FROM TARGET
% ========================================================================

predictedDifferenceFromTarget = ...
    predictedOutwardForPlot - ...
    targetOutwardPositionForPlot;


actualDifferenceFromTarget = ...
    actualOutwardForPlot - ...
    targetOutwardPositionForPlot;


%% Actual-minus-predicted shift
%
% Positive:
%   actual endpoint is farther outward than predicted.
%
% Negative:
%   actual endpoint is less outward than predicted.

actualMinusPredicted = ...
    actualDifferenceFromTarget - ...
    predictedDifferenceFromTarget;


%% ========================================================================
% GRAND REACH-LEVEL MEANS AND SEM
% ========================================================================

grandMeanPredicted = ...
    mean(predictedDifferenceFromTarget,'omitnan');


grandMeanActual = ...
    mean(actualDifferenceFromTarget,'omitnan');


semPredicted = ...
    std(predictedDifferenceFromTarget,0,'omitnan') ./ ...
    sqrt(numberValidReaches);


semActual = ...
    std(actualDifferenceFromTarget,0,'omitnan') ./ ...
    sqrt(numberValidReaches);


grandMeanShift = ...
    mean(actualMinusPredicted,'omitnan');


semShift = ...
    std(actualMinusPredicted,0,'omitnan') ./ ...
    sqrt(numberValidReaches);


%% ========================================================================
% REACH-LEVEL PAIRED T-TEST
%
% This test treats reaches as independent observations.
%
% Because multiple reaches come from the same animal, this p-value should
% be considered descriptive or supplementary. Animal-level statistics or
% a mixed-effects model should be used for primary inference.
% ========================================================================

[hReach,pReach,confidenceIntervalReach,reachTTestStats] = ...
    ttest( ...
        predictedDifferenceFromTarget, ...
        actualDifferenceFromTarget);


%% Optional one-sided test:
%
% This tests whether actual values are greater than predicted values along
% the outward axis.

[~,pActualGreaterThanPredicted] = ...
    ttest( ...
        actualDifferenceFromTarget, ...
        predictedDifferenceFromTarget, ...
        'Tail','right');


%% ========================================================================
% PRINT RESULTS
% ========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('REACH-LEVEL SIGNED DIFFERENCE FROM TARGET\n');
fprintf('ALL CSM SUBTYPES COMBINED\n');
fprintf('============================================================\n');


fprintf('Paired reaches: %d\n', ...
    numberValidReaches);


fprintf('\nPredicted difference from target: %.5f +/- %.5f %s\n', ...
    grandMeanPredicted, ...
    semPredicted, ...
    positionUnit);


fprintf('Actual difference from target: %.5f +/- %.5f %s\n', ...
    grandMeanActual, ...
    semActual, ...
    positionUnit);


fprintf('Actual - predicted shift: %.5f +/- %.5f %s\n', ...
    grandMeanShift, ...
    semShift, ...
    positionUnit);


fprintf('Two-sided reach-level paired t(%d) = %.4f, p = %.6g\n', ...
    reachTTestStats.df, ...
    reachTTestStats.tstat, ...
    pReach);


fprintf('One-sided actual > predicted p = %.6g\n', ...
    pActualGreaterThanPredicted);


fprintf(['95%% CI of predicted - actual difference: ' ...
    '[%.5f, %.5f] %s\n'], ...
    confidenceIntervalReach(1), ...
    confidenceIntervalReach(2), ...
    positionUnit);


%% ========================================================================
% PLOT
% ========================================================================

figureReachSignedDifference = figure( ...
    'Name','Reach-level signed target difference', ...
    'Color','w', ...
    'Position',figurePosition);


axReach = axes(figureReachSignedDifference);

hold(axReach,'on');


%% ------------------------------------------------------------------------
% Draw every paired reach in light gray first
%
% These are drawn before all points so they remain in the background.
% -------------------------------------------------------------------------

for reachIdx = 1:numberValidReaches

    plot( ...
        axReach, ...
        [1 2], ...
        [ ...
            predictedDifferenceFromTarget(reachIdx), ...
            actualDifferenceFromTarget(reachIdx)], ...
        '-', ...
        'Color',pairedLineColor, ...
        'LineWidth',pairedLineWidth, ...
        'HandleVisibility','off');

end


%% ------------------------------------------------------------------------
% Plot every predicted reach in dark gray
% -------------------------------------------------------------------------
predictedPointColor1=rgb('Gray');
predictedReachHandle = scatter( ...
    axReach, ...
    ones(numberValidReaches,1), ...
    predictedDifferenceFromTarget, ...
    predictedPointSize, ...
    predictedPointColor1, ...
    'filled', ...
    'MarkerFaceAlpha',individualPointAlpha, ...
    'MarkerEdgeAlpha',0.10);


%% ------------------------------------------------------------------------
% Plot every actual reach in purple
% -------------------------------------------------------------------------
actualPointColor1=rgb('Gray');
actualReachHandle = scatter( ...
    axReach, ...
    2 .* ones(numberValidReaches,1), ...
    actualDifferenceFromTarget, ...
    actualPointSize, ...
    actualPointColor1, ...
    'filled', ...
    'MarkerFaceAlpha',individualPointAlpha, ...
    'MarkerEdgeAlpha',0.10);


%% ------------------------------------------------------------------------
% Draw grand mean connection
% -------------------------------------------------------------------------

grandMeanConnectionHandle = plot( ...
    axReach, ...
    [1 2], ...
    [ ...
        grandMeanPredicted, ...
        grandMeanActual], ...
    'k-', ...
    'LineWidth',2.3);


%% ------------------------------------------------------------------------
% Predicted grand mean +/- SEM
% -------------------------------------------------------------------------

predictedGrandMeanHandle = errorbar( ...
    axReach, ...
    1, ...
    grandMeanPredicted, ...
    semPredicted, ...
    'o', ...
    'Color',predictedPointColor, ...
    'MarkerFaceColor',predictedPointColor, ...
    'MarkerEdgeColor','k', ...
    'MarkerSize',grandMeanMarkerSize, ...
    'LineWidth',grandMeanLineWidth, ...
    'CapSize',12);


%% ------------------------------------------------------------------------
% Actual grand mean +/- SEM
% -------------------------------------------------------------------------

actualGrandMeanHandle = errorbar( ...
    axReach, ...
    2, ...
    grandMeanActual, ...
    semActual, ...
    'o', ...
    'Color',actualPointColor, ...
    'MarkerFaceColor',actualPointColor, ...
    'MarkerEdgeColor','k', ...
    'MarkerSize',grandMeanMarkerSize, ...
    'LineWidth',grandMeanLineWidth, ...
    'CapSize',12);


%% Redraw grand mean points over the connection line

plot( ...
    axReach, ...
    1, ...
    grandMeanPredicted, ...
    'o', ...
    'MarkerSize',grandMeanMarkerSize, ...
    'MarkerFaceColor',predictedPointColor, ...
    'MarkerEdgeColor','k', ...
    'LineWidth',1.4, ...
    'HandleVisibility','off');


plot( ...
    axReach, ...
    2, ...
    grandMeanActual, ...
    'o', ...
    'MarkerSize',grandMeanMarkerSize, ...
    'MarkerFaceColor',actualPointColor, ...
    'MarkerEdgeColor','k', ...
    'LineWidth',1.4, ...
    'HandleVisibility','off');


%% ------------------------------------------------------------------------
% Target line
%
% Since all values are target-centered, zero represents the target.
% -------------------------------------------------------------------------

targetLineHandle = yline( ...
    axReach, ...
    0, ...
    'k--', ...
    'Target', ...
    'LineWidth',1.7, ...
    'LabelHorizontalAlignment','right', ...
    'LabelVerticalAlignment','bottom');


%% ========================================================================
% AXIS FORMATTING
% ========================================================================

xlim(axReach,[0.68 2.32]);


xticks(axReach,[1 2]);


xticklabels(axReach,{ ...
    'Predicted', ...
    'Actual'});


ylabel(axReach,sprintf( ...
    'Signed outward difference from target (%s)', ...
    positionUnit));


title(axReach,{ ...
    sprintf('All CSM reaches, n = %d', ...
        numberValidReaches), ...
    sprintf('Reach-level paired t(%d) = %.2f, p = %.3g', ...
        reachTTestStats.df, ...
        reachTTestStats.tstat, ...
        pReach)});


box(axReach,'off');


set( ...
    axReach, ...
    'FontSize',12, ...
    'LineWidth',1.3, ...
    'TickDir','out');


%% ========================================================================
% Y-LIMITS
% ========================================================================

if all(isfinite(manualYLimits))

    ylim(axReach,manualYLimits);

else

    maximumAbsoluteDifference = ...
        max(abs([ ...
            predictedDifferenceFromTarget; ...
            actualDifferenceFromTarget]));


    if ~isfinite(maximumAbsoluteDifference) || ...
            maximumAbsoluteDifference == 0

        maximumAbsoluteDifference = 1;

    end


    ylim(axReach,[ ...
        -1.10 .* maximumAbsoluteDifference, ...
         1.10 .* maximumAbsoluteDifference]);

end
ylim([-1.75 1])

%% ========================================================================
% LEGEND
% ========================================================================

% legend( ...
%     axReach, ...
%     [ ...
%         predictedReachHandle, ...
%         actualReachHandle, ...
%         predictedGrandMeanHandle, ...
%         actualGrandMeanHandle, ...
%         targetLineHandle], ...
%     { ...
%         'Individual predicted reaches', ...
%         'Individual actual reaches', ...
%         'Predicted grand mean +/- SEM', ...
%         'Actual grand mean +/- SEM', ...
%         'Target'}, ...
%     'Location','bestoutside');
% 

%% ========================================================================
% CREATE REACH-LEVEL TABLE
% ========================================================================

ReachIndex = ...
    (1:numberValidReaches)';


AnimalID = ...
    allAnimalIDs;


SubtypeID = ...
    allSubtypeIDs;


Subtype = ...
    allSubtypeLabels;


OriginalReachIndex = ...
    allOriginalReachIndices;


PredictedOutwardPosition = ...
    predictedOutwardForPlot;


ActualOutwardPosition = ...
    actualOutwardForPlot;


PredictedDifferenceFromTarget = ...
    predictedDifferenceFromTarget;


ActualDifferenceFromTarget = ...
    actualDifferenceFromTarget;


ActualMinusPredicted = ...
    actualMinusPredicted;


ReachLevelSignedDifferenceTable = table( ...
    ReachIndex, ...
    AnimalID, ...
    SubtypeID, ...
    Subtype, ...
    OriginalReachIndex, ...
    PredictedOutwardPosition, ...
    ActualOutwardPosition, ...
    PredictedDifferenceFromTarget, ...
    ActualDifferenceFromTarget, ...
    ActualMinusPredicted);


disp(ReachLevelSignedDifferenceTable);


%% ========================================================================
% STORE RESULTS
% ========================================================================

ReachLevelSignedDifferenceResults = struct;


%% Settings

ReachLevelSignedDifferenceResults.settings.subtypes = ...
    combinedSubtypes;

ReachLevelSignedDifferenceResults.settings.outwardEndpointColumn = ...
    outwardEndpointColumn;

ReachLevelSignedDifferenceResults.settings.positionScale = ...
    positionScale;

ReachLevelSignedDifferenceResults.settings.positionUnit = ...
    positionUnit;

ReachLevelSignedDifferenceResults.settings.targetOutwardPosition = ...
    targetOutwardPositionForPlot;


%% Reach-level results

ReachLevelSignedDifferenceResults.reach.predictedOutwardPosition = ...
    predictedOutwardForPlot;

ReachLevelSignedDifferenceResults.reach.actualOutwardPosition = ...
    actualOutwardForPlot;

ReachLevelSignedDifferenceResults.reach.predictedDifference = ...
    predictedDifferenceFromTarget;

ReachLevelSignedDifferenceResults.reach.actualDifference = ...
    actualDifferenceFromTarget;

ReachLevelSignedDifferenceResults.reach.actualMinusPredicted = ...
    actualMinusPredicted;

ReachLevelSignedDifferenceResults.reach.animalID = ...
    allAnimalIDs;

ReachLevelSignedDifferenceResults.reach.subtypeID = ...
    allSubtypeIDs;

ReachLevelSignedDifferenceResults.reach.subtypeLabel = ...
    allSubtypeLabels;

ReachLevelSignedDifferenceResults.reach.originalReachIndex = ...
    allOriginalReachIndices;


%% Grand reach-level results

ReachLevelSignedDifferenceResults.group.predictedMean = ...
    grandMeanPredicted;

ReachLevelSignedDifferenceResults.group.actualMean = ...
    grandMeanActual;

ReachLevelSignedDifferenceResults.group.predictedSEM = ...
    semPredicted;

ReachLevelSignedDifferenceResults.group.actualSEM = ...
    semActual;

ReachLevelSignedDifferenceResults.group.meanActualMinusPredicted = ...
    grandMeanShift;

ReachLevelSignedDifferenceResults.group.semActualMinusPredicted = ...
    semShift;


%% Statistics

ReachLevelSignedDifferenceResults.statistics.h = ...
    hReach;

ReachLevelSignedDifferenceResults.statistics.pTwoSided = ...
    pReach;

ReachLevelSignedDifferenceResults.statistics. ...
    pActualGreaterThanPredicted = ...
    pActualGreaterThanPredicted;

ReachLevelSignedDifferenceResults.statistics.t = ...
    reachTTestStats.tstat;

ReachLevelSignedDifferenceResults.statistics.df = ...
    reachTTestStats.df;

ReachLevelSignedDifferenceResults.statistics.confidenceInterval = ...
    confidenceIntervalReach;


%% Table

ReachLevelSignedDifferenceResults.reachTable = ...
    ReachLevelSignedDifferenceTable;


fprintf('\nAnalysis complete.\n');
fprintf('Results stored in ReachLevelSignedDifferenceResults.\n');

%% ========================================================================
% PREDICTED VERSUS ACTUAL SIGNED DIFFERENCE FROM TARGET
% ALL CSM SUBTYPES COMBINED
%
% VISUALIZATION
%
%   Light-gray lines:
%       Individual predicted-to-actual reach pairs
%
%   Dark-gray points:
%       Individual predicted endpoints
%
%   Purple points:
%       Individual actual endpoints
%
%   Large dark-gray and purple points:
%       Grand means calculated across animal means
%
%   Error bars:
%       SEM calculated across animal means
%
%   Dashed line at zero:
%       Target position
%
% SIGNED DIFFERENCE FROM TARGET
%
%   endpoint outward position - target outward position
%
% Interpretation:
%
%   0        = endpoint is at the target
%   negative = endpoint falls short of the target
%   positive = endpoint extends beyond the target
%
% PRIMARY STATISTICS
%
% A paired t-test is performed on one predicted and one actual mean per
% animal. This is the appropriate inferential test because reaches are
% nested within animals.
%
% REQUIRED WORKSPACE VARIABLES
%
%   DecelerationEndpointResults
%   targetPosition
%   analysisDimensions
%
% Expected subtype fields:
%
%   DecelerationEndpointResults.T0.reach.predictedEndpoint
%   DecelerationEndpointResults.T0.reach.actualEndpoint
%   DecelerationEndpointResults.T0.reach.animalID
%
% and corresponding fields for T1, T2, and T3.
% ========================================================================


%% ========================================================================
% SETTINGS
% ========================================================================

combinedSubtypes = {'T0','T1','T2','T3'};


%% Endpoint coordinate to analyze
%
% With:
%
%   analysisDimensions = [1 2]
%
% endpoint column 1 = x
% endpoint column 2 = y
%
% Use 1 for outward x position.

outwardEndpointColumn = 1;


%% Units
%
% Original positions are in decimeters.
%
% For centimeters:

positionScale = 10;
positionUnit = 'cm';

% For decimeters:
%
% positionScale = 1;
% positionUnit = 'dm';


%% Animal inclusion threshold
%
% Each animal must contribute at least this many valid reaches after T0-T3
% are combined.

minimumReachesPerAnimal = 5;


%% Plot colors

predictedPointColor = [0.25 0.25 0.25];

actualPointColor = [0.5781 0 0.8242];

pairedLineColor = [0.82 0.82 0.82];


%% Plot sizes

individualPointSize = 18;

individualPointAlpha = 0.32;

pairedLineWidth = 0.40;

grandMeanPointSize = 145;

grandMeanErrorLineWidth = 3.0;

grandMeanCapSize = 18;


%% Wider figure

%figurePosition = [180 100 650 700];
%figurePosition = [200 200 420 900];
figurePosition = [200 200 400 350];


%% Optional manual y-limits
%
% Leave as NaN for automatic symmetric limits around zero.
%
% Example:
%
% manualYLimits = [-1.5 1.5];

manualYLimits = [NaN NaN];


%% ========================================================================
% VALIDATE INPUTS
% ========================================================================

if ~exist('DecelerationEndpointResults','var')

    error('DecelerationEndpointResults is missing from the workspace.');

end


if ~exist('targetPosition','var') || ...
        isempty(targetPosition)

    error('targetPosition is missing from the workspace.');

end


if ~exist('analysisDimensions','var') || ...
        isempty(analysisDimensions)

    error('analysisDimensions is missing from the workspace.');

end


if outwardEndpointColumn > numel(analysisDimensions)

    error(['outwardEndpointColumn exceeds the number of coordinates ' ...
        'stored in the endpoint arrays.']);

end


%% ========================================================================
% TARGET OUTWARD POSITION
% ========================================================================

targetCoordinatesUsed = ...
    targetPosition(analysisDimensions);


targetOutwardOriginalUnits = ...
    targetCoordinatesUsed(outwardEndpointColumn);


targetOutwardForPlot = ...
    targetOutwardOriginalUnits .* positionScale;


fprintf('\nTarget outward coordinate: %.5f %s\n', ...
    targetOutwardForPlot, ...
    positionUnit);


%% ========================================================================
% COMBINE VALID T0-T3 REACHES
% ========================================================================

allPredictedOutward = [];

allActualOutward = [];

allAnimalIDs = [];

allSubtypeIDs = [];

allSubtypeLabels = {};

allOriginalReachIndices = [];


for subtypeIdx = 1:numel(combinedSubtypes)

    subtype = combinedSubtypes{subtypeIdx};


    if ~isfield(DecelerationEndpointResults,subtype)

        warning('%s is missing from DecelerationEndpointResults.', ...
            subtype);

        continue

    end


    R = DecelerationEndpointResults.(subtype);


    if ~isfield(R,'reach') || ...
            ~isfield(R.reach,'predictedEndpoint') || ...
            ~isfield(R.reach,'actualEndpoint') || ...
            ~isfield(R.reach,'animalID')

        warning('%s does not contain the required reach-level fields.', ...
            subtype);

        continue

    end


    predictedEndpoints = ...
        R.reach.predictedEndpoint;

    actualEndpoints = ...
        R.reach.actualEndpoint;

    animalIDs = ...
        R.reach.animalID;


    if isfield(R.reach,'originalReachIndex')

        originalReachIndices = ...
            R.reach.originalReachIndex;

    else

        originalReachIndices = ...
            (1:size(predictedEndpoints,1))';

    end


    %% Ensure all arrays have matching lengths

    numberRows = min([ ...
        size(predictedEndpoints,1), ...
        size(actualEndpoints,1), ...
        numel(animalIDs), ...
        numel(originalReachIndices)]);


    predictedEndpoints = ...
        predictedEndpoints(1:numberRows,:);

    actualEndpoints = ...
        actualEndpoints(1:numberRows,:);

    animalIDs = ...
        animalIDs(1:numberRows);

    originalReachIndices = ...
        originalReachIndices(1:numberRows);


    %% Keep valid matched endpoint pairs

    validPairs = ...
        all(isfinite(predictedEndpoints),2) & ...
        all(isfinite(actualEndpoints),2) & ...
        isfinite(animalIDs);


    predictedEndpoints = ...
        predictedEndpoints(validPairs,:);

    actualEndpoints = ...
        actualEndpoints(validPairs,:);

    animalIDs = ...
        animalIDs(validPairs);

    originalReachIndices = ...
        originalReachIndices(validPairs);


    %% Extract outward coordinate

    predictedOutward = ...
        predictedEndpoints(:,outwardEndpointColumn);

    actualOutward = ...
        actualEndpoints(:,outwardEndpointColumn);


    validOutward = ...
        isfinite(predictedOutward) & ...
        isfinite(actualOutward);


    predictedOutward = ...
        predictedOutward(validOutward);

    actualOutward = ...
        actualOutward(validOutward);

    animalIDs = ...
        animalIDs(validOutward);

    originalReachIndices = ...
        originalReachIndices(validOutward);


    numberSubtypeReaches = ...
        numel(predictedOutward);


    %% Add subtype data to combined arrays

    allPredictedOutward = [ ...
        allPredictedOutward; ...
        predictedOutward];


    allActualOutward = [ ...
        allActualOutward; ...
        actualOutward];


    allAnimalIDs = [ ...
        allAnimalIDs; ...
        animalIDs];


    allSubtypeIDs = [ ...
        allSubtypeIDs; ...
        repmat(subtypeIdx,numberSubtypeReaches,1)];


    allSubtypeLabels = [ ...
        allSubtypeLabels; ...
        repmat({subtype},numberSubtypeReaches,1)];


    allOriginalReachIndices = [ ...
        allOriginalReachIndices; ...
        originalReachIndices];


    fprintf('%s: added %d valid reaches\n', ...
        subtype, ...
        numberSubtypeReaches);

end


numberValidReaches = ...
    numel(allPredictedOutward);


if numberValidReaches < 2

    error('Fewer than two valid predicted-actual reach pairs were found.');

end


fprintf('\nTotal valid paired reaches: %d\n', ...
    numberValidReaches);


%% ========================================================================
% CONVERT TO DISPLAY UNITS
% ========================================================================

predictedOutwardForPlot = ...
    allPredictedOutward .* positionScale;


actualOutwardForPlot = ...
    allActualOutward .* positionScale;


%% ========================================================================
% SIGNED DIFFERENCE FROM TARGET FOR EVERY REACH
% ========================================================================

predictedDifferenceFromTarget = ...
    predictedOutwardForPlot - targetOutwardForPlot;


actualDifferenceFromTarget = ...
    actualOutwardForPlot - targetOutwardForPlot;


actualMinusPredicted = ...
    actualDifferenceFromTarget - ...
    predictedDifferenceFromTarget;


%% ========================================================================
% ANIMAL-LEVEL MEANS
% ========================================================================

uniqueAnimalIDs = ...
    unique(allAnimalIDs);


numberCandidateAnimals = ...
    numel(uniqueAnimalIDs);


animalPredictedMean = ...
    nan(numberCandidateAnimals,1);


animalActualMean = ...
    nan(numberCandidateAnimals,1);


animalReachCount = ...
    zeros(numberCandidateAnimals,1);


animalSubtypeReachCount = ...
    zeros(numberCandidateAnimals,numel(combinedSubtypes));


for animalIdx = 1:numberCandidateAnimals

    thisAnimalID = ...
        uniqueAnimalIDs(animalIdx);


    animalMask = ...
        allAnimalIDs == thisAnimalID;


    animalReachCount(animalIdx) = ...
        sum(animalMask);


    for subtypeIdx = 1:numel(combinedSubtypes)

        animalSubtypeReachCount(animalIdx,subtypeIdx) = ...
            sum( ...
                animalMask & ...
                allSubtypeIDs == subtypeIdx);

    end


    if animalReachCount(animalIdx) < ...
            minimumReachesPerAnimal

        continue

    end


    animalPredictedMean(animalIdx) = ...
        mean( ...
            predictedDifferenceFromTarget(animalMask), ...
            'omitnan');


    animalActualMean(animalIdx) = ...
        mean( ...
            actualDifferenceFromTarget(animalMask), ...
            'omitnan');

end


%% Valid paired animals

validAnimals = ...
    isfinite(animalPredictedMean) & ...
    isfinite(animalActualMean);


pairedAnimalIDs = ...
    uniqueAnimalIDs(validAnimals);


pairedAnimalPredicted = ...
    animalPredictedMean(validAnimals);


pairedAnimalActual = ...
    animalActualMean(validAnimals);


pairedAnimalReachCount = ...
    animalReachCount(validAnimals);


pairedAnimalSubtypeReachCount = ...
    animalSubtypeReachCount(validAnimals,:);


numberValidAnimals = ...
    sum(validAnimals);


if numberValidAnimals < 2

    error(['Fewer than two animals contributed sufficient paired ' ...
        'predicted and actual reaches.']);

end


%% ========================================================================
% GRAND MEAN AND SEM ACROSS ANIMALS
% ========================================================================

grandMeanPredicted = ...
    mean(pairedAnimalPredicted,'omitnan');


grandMeanActual = ...
    mean(pairedAnimalActual,'omitnan');


semPredicted = ...
    std(pairedAnimalPredicted,0,'omitnan') ./ ...
    sqrt(numberValidAnimals);


semActual = ...
    std(pairedAnimalActual,0,'omitnan') ./ ...
    sqrt(numberValidAnimals);


animalActualMinusPredicted = ...
    pairedAnimalActual - pairedAnimalPredicted;


grandMeanShift = ...
    mean(animalActualMinusPredicted,'omitnan');


semShift = ...
    std(animalActualMinusPredicted,0,'omitnan') ./ ...
    sqrt(numberValidAnimals);


%% ========================================================================
% PRIMARY ANIMAL-LEVEL PAIRED T-TEST
% ========================================================================

[hAnimal,pAnimal,ciAnimal,animalTTestStats] = ...
    ttest( ...
        pairedAnimalPredicted, ...
        pairedAnimalActual);


%% Optional directional test
%
% Tests whether actual signed outward position is greater than predicted.

[~,pAnimalActualGreater] = ...
    ttest( ...
        pairedAnimalActual, ...
        pairedAnimalPredicted, ...
        'Tail','right');


%% ========================================================================
% DESCRIPTIVE REACH-LEVEL PAIRED T-TEST
%
% This treats reaches as independent and should not be the primary
% inferential statistic.
% ========================================================================

[hReach,pReach,ciReach,reachTTestStats] = ...
    ttest( ...
        predictedDifferenceFromTarget, ...
        actualDifferenceFromTarget);


%% ========================================================================
% PRINT RESULTS
% ========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('SIGNED DIFFERENCE FROM TARGET\n');
fprintf('ALL CSM SUBTYPES COMBINED\n');
fprintf('============================================================\n');


fprintf('Valid animals: %d\n', ...
    numberValidAnimals);


fprintf('Valid reaches: %d\n', ...
    numberValidReaches);


fprintf('\nAnimal-level predicted mean: %.5f +/- %.5f %s\n', ...
    grandMeanPredicted, ...
    semPredicted, ...
    positionUnit);


fprintf('Animal-level actual mean: %.5f +/- %.5f %s\n', ...
    grandMeanActual, ...
    semActual, ...
    positionUnit);


fprintf('Animal-level actual - predicted shift: %.5f +/- %.5f %s\n', ...
    grandMeanShift, ...
    semShift, ...
    positionUnit);


fprintf('\nPrimary paired animal-level t-test:\n');

fprintf('t(%d) = %.4f, p = %.6g\n', ...
    animalTTestStats.df, ...
    animalTTestStats.tstat, ...
    pAnimal);


fprintf('One-sided actual > predicted p = %.6g\n', ...
    pAnimalActualGreater);


fprintf(['95%% CI of predicted - actual animal-level difference: ' ...
    '[%.5f, %.5f] %s\n'], ...
    ciAnimal(1), ...
    ciAnimal(2), ...
    positionUnit);


fprintf('\nDescriptive reach-level paired test:\n');

fprintf('t(%d) = %.4f, p = %.6g\n', ...
    reachTTestStats.df, ...
    reachTTestStats.tstat, ...
    pReach);


fprintf('\nSEM values plotted across animals:\n');

fprintf('Predicted SEM = %.8f %s\n', ...
    semPredicted, ...
    positionUnit);

fprintf('Actual SEM = %.8f %s\n', ...
    semActual, ...
    positionUnit);


%% ========================================================================
% PLOT
% ========================================================================

figureReachSignedDifference = figure( ...
    'Name','Reach-level signed target difference', ...
    'Color','w', ...
    'Position',figurePosition);


axReach = axes(figureReachSignedDifference);

hold(axReach,'on');


%% ------------------------------------------------------------------------
% Draw individual paired reaches first so they remain in the background
% -------------------------------------------------------------------------

for reachIdx = 1:numberValidReaches

    plot( ...
        axReach, ...
        [1 2], ...
        [ ...
            predictedDifferenceFromTarget(reachIdx), ...
            actualDifferenceFromTarget(reachIdx)], ...
        '-', ...
        'Color',pairedLineColor, ...
        'LineWidth',pairedLineWidth, ...
        'HandleVisibility','off');

end


%% ------------------------------------------------------------------------
% Individual predicted reaches
% -------------------------------------------------------------------------
predictedPointColor1=rgb('gray');
predictedReachHandle = scatter( ...
    axReach, ...
    ones(numberValidReaches,1), ...
    predictedDifferenceFromTarget, ...
    individualPointSize, ...
    predictedPointColor1, ...
    'filled', ...
    'MarkerFaceAlpha',individualPointAlpha, ...
    'MarkerEdgeAlpha',0.10);


%% ------------------------------------------------------------------------
% Individual actual reaches
% -------------------------------------------------------------------------
actualPointColor1=rgb('Gray');
actualReachHandle = scatter( ...
    axReach, ...
    2 .* ones(numberValidReaches,1), ...
    actualDifferenceFromTarget, ...
    individualPointSize, ...
    actualPointColor, ...
    'filled', ...
    'MarkerFaceAlpha',individualPointAlpha, ...
    'MarkerEdgeAlpha',0.10);


%% ------------------------------------------------------------------------
% Draw grand-mean connection first
% -------------------------------------------------------------------------

grandMeanConnectionHandle = plot( ...
    axReach, ...
    [1 2], ...
    [ ...
        grandMeanPredicted, ...
        grandMeanActual], ...
    'k-', ...
    'LineWidth',2.4);


%% ------------------------------------------------------------------------
% Draw SEM bars without markers
%
% This ensures the mean markers do not cover the vertical error bars.
% -------------------------------------------------------------------------

predictedSEMHandle = errorbar( ...
    axReach, ...
    1, ...
    grandMeanPredicted, ...
    semPredicted, ...
    semPredicted, ...
    'LineStyle','none', ...
    'Color',predictedPointColor, ...
    'LineWidth',grandMeanErrorLineWidth, ...
    'CapSize',grandMeanCapSize);


actualSEMHandle = errorbar( ...
    axReach, ...
    2, ...
    grandMeanActual, ...
    semActual, ...
    semActual, ...
    'LineStyle','none', ...
    'Color',actualPointColor, ...
    'LineWidth',grandMeanErrorLineWidth, ...
    'CapSize',grandMeanCapSize);


%% ------------------------------------------------------------------------
% Draw grand-mean points last so they appear on top
% -------------------------------------------------------------------------

predictedGrandMeanHandle = scatter( ...
    axReach, ...
    1, ...
    grandMeanPredicted, ...
    grandMeanPointSize, ...
    predictedPointColor, ...
    'filled', ...
    'MarkerEdgeColor','k', ...
    'LineWidth',1.5);


actualGrandMeanHandle = scatter( ...
    axReach, ...
    2, ...
    grandMeanActual, ...
    grandMeanPointSize, ...
    actualPointColor, ...
    'filled', ...
    'MarkerEdgeColor','k', ...
    'LineWidth',1.5);


%% ------------------------------------------------------------------------
% Target line
% -------------------------------------------------------------------------

targetLineHandle = yline( ...
    axReach, ...
    0, ...
    'k--', ...
    'Target', ...
    'LineWidth',1.7, ...
    'LabelHorizontalAlignment','right', ...
    'LabelVerticalAlignment','bottom');


%% ========================================================================
% FORMATTING
% ========================================================================

xlim(axReach,[0.68 2.32]);


xticks(axReach,[1 2]);


xticklabels(axReach,{ ...
    'Predicted', ...
    'Actual'});


ylabel(axReach,sprintf( ...
    'Signed outward difference from target (%s)', ...
    positionUnit));


title(axReach,{ ...
    sprintf('All CSM reaches: n = %d reaches, N = %d animals', ...
        numberValidReaches, ...
        numberValidAnimals), ...
    sprintf('Animal-level paired t(%d) = %.2f, p = %.3g', ...
        animalTTestStats.df, ...
        animalTTestStats.tstat, ...
        pAnimal)});


box(axReach,'off');


set( ...
    axReach, ...
    'FontSize',12, ...
    'LineWidth',1.3, ...
    'TickDir','out');


%% ========================================================================
% Y-LIMITS
% ========================================================================

if all(isfinite(manualYLimits))

    ylim(axReach,manualYLimits);

else

    allPlotValues = [ ...
        predictedDifferenceFromTarget; ...
        actualDifferenceFromTarget; ...
        grandMeanPredicted-semPredicted; ...
        grandMeanPredicted+semPredicted; ...
        grandMeanActual-semActual; ...
        grandMeanActual+semActual; ...
        0];


    maximumAbsoluteDifference = ...
        max(abs(allPlotValues));


    if ~isfinite(maximumAbsoluteDifference) || ...
            maximumAbsoluteDifference == 0

        maximumAbsoluteDifference = 1;

    end


    ylim(axReach,[ ...
        -1.10 .* maximumAbsoluteDifference, ...
         1.10 .* maximumAbsoluteDifference]);

end
ylim([-1.75 1])

%% ========================================================================
% LEGEND
% ========================================================================

% legend( ...
%     axReach, ...
%     [ ...
%         predictedReachHandle, ...
%         actualReachHandle, ...
%         predictedGrandMeanHandle, ...
%         actualGrandMeanHandle, ...
%         targetLineHandle], ...
%     { ...
%         'Individual predicted reaches', ...
%         'Individual actual reaches', ...
%         'Predicted grand mean +/- SEM across animals', ...
%         'Actual grand mean +/- SEM across animals', ...
%         'Target'}, ...
%     'Location','bestoutside');


%% ========================================================================
% REACH-LEVEL TABLE
% ========================================================================

ReachIndex = ...
    (1:numberValidReaches)';


AnimalID = ...
    allAnimalIDs;


SubtypeID = ...
    allSubtypeIDs;


Subtype = ...
    allSubtypeLabels;


OriginalReachIndex = ...
    allOriginalReachIndices;


PredictedOutwardPosition = ...
    predictedOutwardForPlot;


ActualOutwardPosition = ...
    actualOutwardForPlot;


PredictedDifferenceFromTarget = ...
    predictedDifferenceFromTarget;


ActualDifferenceFromTarget = ...
    actualDifferenceFromTarget;


ActualMinusPredicted = ...
    actualMinusPredicted;


ReachLevelSignedDifferenceTable = table( ...
    ReachIndex, ...
    AnimalID, ...
    SubtypeID, ...
    Subtype, ...
    OriginalReachIndex, ...
    PredictedOutwardPosition, ...
    ActualOutwardPosition, ...
    PredictedDifferenceFromTarget, ...
    ActualDifferenceFromTarget, ...
    ActualMinusPredicted);


disp(ReachLevelSignedDifferenceTable);


%% ========================================================================
% ANIMAL-LEVEL TABLE
% ========================================================================

if exist('animalNames','var') && ...
        numel(animalNames) >= max(pairedAnimalIDs)

    Animal = ...
        animalNames(pairedAnimalIDs);

else

    Animal = ...
        compose('Animal_%d',pairedAnimalIDs);

end


T0ReachCount = ...
    pairedAnimalSubtypeReachCount(:,1);

T1ReachCount = ...
    pairedAnimalSubtypeReachCount(:,2);

T2ReachCount = ...
    pairedAnimalSubtypeReachCount(:,3);

T3ReachCount = ...
    pairedAnimalSubtypeReachCount(:,4);

TotalReachCount = ...
    pairedAnimalReachCount;

PredictedAnimalMean = ...
    pairedAnimalPredicted;

ActualAnimalMean = ...
    pairedAnimalActual;

ActualMinusPredictedAnimalMean = ...
    animalActualMinusPredicted;


AnimalSignedDifferenceTable = table( ...
    Animal, ...
    T0ReachCount, ...
    T1ReachCount, ...
    T2ReachCount, ...
    T3ReachCount, ...
    TotalReachCount, ...
    PredictedAnimalMean, ...
    ActualAnimalMean, ...
    ActualMinusPredictedAnimalMean);


disp(AnimalSignedDifferenceTable);


%% ========================================================================
% STORE RESULTS
% ========================================================================

ReachLevelSignedDifferenceResults = struct;


%% Settings

ReachLevelSignedDifferenceResults.settings.subtypes = ...
    combinedSubtypes;

ReachLevelSignedDifferenceResults.settings.outwardEndpointColumn = ...
    outwardEndpointColumn;

ReachLevelSignedDifferenceResults.settings.positionScale = ...
    positionScale;

ReachLevelSignedDifferenceResults.settings.positionUnit = ...
    positionUnit;

ReachLevelSignedDifferenceResults.settings.targetOutwardPosition = ...
    targetOutwardForPlot;

ReachLevelSignedDifferenceResults.settings.minimumReachesPerAnimal = ...
    minimumReachesPerAnimal;


%% Reach-level data

ReachLevelSignedDifferenceResults.reach.predictedOutwardPosition = ...
    predictedOutwardForPlot;

ReachLevelSignedDifferenceResults.reach.actualOutwardPosition = ...
    actualOutwardForPlot;

ReachLevelSignedDifferenceResults.reach.predictedDifference = ...
    predictedDifferenceFromTarget;

ReachLevelSignedDifferenceResults.reach.actualDifference = ...
    actualDifferenceFromTarget;

ReachLevelSignedDifferenceResults.reach.actualMinusPredicted = ...
    actualMinusPredicted;

ReachLevelSignedDifferenceResults.reach.animalID = ...
    allAnimalIDs;

ReachLevelSignedDifferenceResults.reach.subtypeID = ...
    allSubtypeIDs;

ReachLevelSignedDifferenceResults.reach.subtypeLabel = ...
    allSubtypeLabels;

ReachLevelSignedDifferenceResults.reach.originalReachIndex = ...
    allOriginalReachIndices;


%% Animal-level data

ReachLevelSignedDifferenceResults.animal.ID = ...
    pairedAnimalIDs;

ReachLevelSignedDifferenceResults.animal.predictedMean = ...
    pairedAnimalPredicted;

ReachLevelSignedDifferenceResults.animal.actualMean = ...
    pairedAnimalActual;

ReachLevelSignedDifferenceResults.animal.actualMinusPredicted = ...
    animalActualMinusPredicted;

ReachLevelSignedDifferenceResults.animal.reachCount = ...
    pairedAnimalReachCount;

ReachLevelSignedDifferenceResults.animal.subtypeReachCount = ...
    pairedAnimalSubtypeReachCount;


%% Group summary

ReachLevelSignedDifferenceResults.group.predictedMean = ...
    grandMeanPredicted;

ReachLevelSignedDifferenceResults.group.actualMean = ...
    grandMeanActual;

ReachLevelSignedDifferenceResults.group.predictedSEM = ...
    semPredicted;

ReachLevelSignedDifferenceResults.group.actualSEM = ...
    semActual;

ReachLevelSignedDifferenceResults.group.meanActualMinusPredicted = ...
    grandMeanShift;

ReachLevelSignedDifferenceResults.group.semActualMinusPredicted = ...
    semShift;


%% Primary animal-level statistics

ReachLevelSignedDifferenceResults.statistics.animal.h = ...
    hAnimal;

ReachLevelSignedDifferenceResults.statistics.animal.pTwoSided = ...
    pAnimal;

ReachLevelSignedDifferenceResults.statistics.animal. ...
    pActualGreaterThanPredicted = ...
    pAnimalActualGreater;

ReachLevelSignedDifferenceResults.statistics.animal.t = ...
    animalTTestStats.tstat;

ReachLevelSignedDifferenceResults.statistics.animal.df = ...
    animalTTestStats.df;

ReachLevelSignedDifferenceResults.statistics.animal. ...
    confidenceInterval = ...
    ciAnimal;


%% Descriptive reach-level statistics

ReachLevelSignedDifferenceResults.statistics.reach.h = ...
    hReach;

ReachLevelSignedDifferenceResults.statistics.reach.p = ...
    pReach;

ReachLevelSignedDifferenceResults.statistics.reach.t = ...
    reachTTestStats.tstat;

ReachLevelSignedDifferenceResults.statistics.reach.df = ...
    reachTTestStats.df;

ReachLevelSignedDifferenceResults.statistics.reach. ...
    confidenceInterval = ...
    ciReach;


%% Tables

ReachLevelSignedDifferenceResults.reachTable = ...
    ReachLevelSignedDifferenceTable;

ReachLevelSignedDifferenceResults.animalTable = ...
    AnimalSignedDifferenceTable;


fprintf('\nAnalysis complete.\n');
fprintf('Results stored in ReachLevelSignedDifferenceResults.\n');