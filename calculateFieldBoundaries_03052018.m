function fieldLocations = calculateFieldBoundaries_03052018(centerOfMass, fieldSize)
% This function calculates the boundaries of the place fields based on the
% average center of mass and the average field width calculated in 
% 'PFAnalysis_LAST_RIGHT_freqnorm_debug.m' 

trackLength = 297;

for i = 1:length(centerOfMass)
    % Correct for center of mass averages greater than the defined track
    % length. The PF analysis function does not properly correct for it.
    if centerOfMass(i) > trackLength
        centerOfMass(i) = centerOfMass(i) - trackLength;
    end
    % Calculate the boundaries of each place field
    beforeCenter = centerOfMass(i) - fieldSize(i)/2;
    afterCenter = centerOfMass(i) + fieldSize(i)/2;
    if beforeCenter < 0
        fieldLocations{i,1} = [0 afterCenter; (trackLength + beforeCenter) trackLength];
    elseif afterCenter > trackLength
        fieldLocations{i,1} = [beforeCenter trackLength; 0 (afterCenter - trackLength)];
    else
        fieldLocations{i,1} = [beforeCenter afterCenter];
    end
end