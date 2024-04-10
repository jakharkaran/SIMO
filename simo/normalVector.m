function [ normalVec ] = normalVector( outNew, normalVec, N, volumeCurr, delVolumeFirst, pointIndex  )
%% normalVector determines appropriate direction of normal vector for
%   perturbation by PN and -PN operators.
%   Outwards normal vector for +PN and inwards for -PN operator

%   Code perturbs a point in the pre-defined direction of normal and
%   calculates volume. If volume increases the normal vector has appropriate
%   direction for +PN operator. If volume decrease the direction of normal
%   vector is reversed. Similarly for -PN operator.

    % Perturbation factor
    P = .1*abs((mean(outNew(1:N/2,1))));

    if pointIndex <= 0
        % Find direction of normal for all points on curve
        % Don't keep pointIndex <=0 unless required to. 

        for countNorm = 1:N/2

            out = outNew;

            % Perturbation vector using pre-defined normal vector
            perturbVec = normalVec(countNorm, 1:2);
            perturbVecSync = [-normalVec(countNorm,1),normalVec(countNorm,2)];
            % Normalise perturbation vector
            perturbVec = perturbVec./norm(perturbVec);
            perturbVecSync = perturbVecSync/norm(perturbVecSync);

            % Perturb point in direction of normal vector
            out(countNorm,1:2) = out(countNorm,1:2) + P*perturbVec; 
            out(N-countNorm+1,1:2) = out(N-countNorm+1,1:2) + P*perturbVecSync;
            out = mvsplint(out,N);
            volumeUpdated = pappus(out);

            % For +PN operator
            if delVolumeFirst > 0
                % Reverse direction of normal vector if volume decreases
                if volumeUpdated < volumeCurr
                    normalVec(countNorm,:) = -normalVec(countNorm,:);
                end
            else
            % For -PN operator
                % Reverse direction of normal vector if volume increases
                if volumeUpdated > volumeCurr
                    normalVec(countNorm,:) = -normalVec(countNorm,:);
                end 
            end
            normalVec(N+1-countNorm,1:2) = [-normalVec(countNorm,1), normalVec(countNorm,2)];
        end
    else
            % Find direction of normal vector on pointIndex

            out = outNew;

            % Perturbation vector using pre-defined normal vector
            perturbVec = normalVec(pointIndex, 1:2);
            perturbVecSync = [-normalVec(pointIndex,1),normalVec(pointIndex,2)];
            % Normalise perturbation vector
            perturbVec = perturbVec./norm(perturbVec);
            perturbVecSync = perturbVecSync/norm(perturbVecSync);

            % Perturb point in direction of normal vector
            out(pointIndex,1:2) = out(pointIndex,1:2) + P*perturbVec; 
            out(N-pointIndex+1,1:2) = out(N-pointIndex+1,1:2) + P*perturbVecSync;
            out = mvsplint(out,N);
            volumeUpdated = pappus(out);

            % For +PN operator
            if delVolumeFirst > 0
                % Reverse direction of normal vector if volume decreases
                if volumeUpdated < volumeCurr
                    normalVec(pointIndex,:) = -normalVec(pointIndex,:);
                end
            % For -PN operator
                % Reverse direction of normal vector if volume increases
            else
                if volumeUpdated > volumeCurr
                    normalVec(pointIndex,:) = -normalVec(pointIndex,:);
                end 
            end
            normalVec(N+1-pointIndex,1:2) = [-normalVec(pointIndex,1), normalVec(pointIndex,2)];
    end
end

