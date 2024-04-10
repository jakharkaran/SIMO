function [outNew, movingInd1, movingInd2, maxIndex] = pin_droplet_sync_add(outNew, N, volumeOrig, normalVec, gravitational, centrifugal)
%%  +PN operator
%    Increases the volume by moving points of maximum curvature in the
%    direction of normal vector, increasing the volume
%    Overall volume change remains constant with error of 10^-5 of the initial volume

    % plt = 1 if you want to see the plots as the code is running 
    plt = 0;

    out = outNew;
    
    % Calculate the current volume
    volumeCurr = pappus(out);

    % Calculate change in volume
    delVolume = volumeOrig - volumeCurr;
    delVolumeFirst = delVolume;

    % Define PID Parameters
    % Parameters are tunable. Following set of parameters deemed suitable for
    % current code. Although of re-tuning the parameters for possibly faster
    % and stable algorithm is not denied.

    P = 50*sqrt(abs(delVolume/2));
    I = P/75;
    D = P/100;

    % Initialize
    errInt = 0;
    errPrev = 0;
    errTol = volumeOrig*1e-5;

    movingInd1 = 1;
    movingInd2 = 1;
    maxInd = 1;
    
    % Initialise no. of loops
    iii = 1;
    %% Compensate volume addition 
    while abs(delVolume) > errTol

        % Index of point with maximum curvature
        [maximum, maxInd] = max(out(2:N-1,4) - out(2:N-1,2)*gravitational + centrifugal*(out(2:N-1,1).^2));
        % End points are constrained
        maxInd = maxInd + 1;

        % Calculating appropriate direction of normal vector of point with
        % maximum curvature for perturbation
        [ normalVec ] = normalVector( out, normalVec, N, volumeCurr, delVolumeFirst, maxInd );
        normalVecUpdated = normalVec;

        % For axisymmetric curve two points symmetric to axis of symmetry are perturbed
        % Perturbation vector
        perturbVec = normalVec(maxInd, 1:2);
        perturbVecSync = normalVec(N+1-maxInd, 1:2);
        % Normalise perturbation vector
        perturbVec = perturbVec./norm(perturbVec);    
        perturbVecSync = perturbVecSync/norm(perturbVecSync);
        
        % Perturb points with minimum curvature to decrease volume
        out(maxInd,1:2) = out(maxInd,1:2) + (P*(delVolume/2)*perturbVec...
            +I*(errInt/2)*perturbVec...
            +D*(delVolume - errPrev)/2*perturbVec);

        out(N-maxInd+1,1:2) = out(N-maxInd+1,1:2) + (P*(delVolume)/2*perturbVecSync...
            +I*(errInt)/2*perturbVecSync...
            +D*(delVolume - errPrev)/2*perturbVecSync);

        % Storing previous volume difference
        errPrev = delVolume;
        
        % Update points with a new spline
        [outUpdated,normalVec] = mvsplint(out,N);    
        out = outUpdated;
        % Calculate the current volume
        volumeCurr=pappus(outUpdated);
        
        % Calculate the new volume difference
        delVolume = volumeOrig - volumeCurr;
        delVolume_arr(iii) = delVolume;
        
        % Compute the error integral term
        errInt = errInt + (delVolume);
        % Initialize the integral term if the error is zero
        if abs(errInt)>abs(volumeCurr)
            errInt = 0;
        end


        %% Output

        outNew = out;

        % Index of perturbing point 
        movingInd1 = maxInd;
        movingInd2 = N-maxInd+1;
        maxIndex(iii) = maxInd;
        %% Plot    
        if plt == 1
            if iii>10
                if mod(iii,20)==0    
                    h = figure(1);
                    subplot(8,2,[1 10]);
                    set(h, 'Position', [50 50 1024 640], 'Color', 'white')
                    plot(out(:,1), out(:,2), 'LineWidth',1.5);
                    axis([-2.5 2.5 0 2]);
                    drawnow
                    daspect([1,1,1]);
                    hold on;

                    quiver(out(maxInd,1),out(maxInd,2),0.1*normalVecUpdated(maxInd,1), .1*normalVecUpdated(maxInd,2));
                    quiver(out(N+1-maxInd,1),out(N+1-maxInd,2),0.1*normalVecUpdated(N+1-maxInd,1), .1*normalVecUpdated(N+1-maxInd,2));

                    delVolume;
                    movingInd1 = maxInd;
                    movingInd2 = N-maxInd+1;

                    plot(out(movingInd1,1), out(movingInd1,2), '-s', 'MarkerEdgeColor','blue','MarkerFaceColor','blue');
                    drawnow           
                    plot(out(movingInd2,1), out(movingInd2,2), '-s', 'MarkerEdgeColor','blue','MarkerFaceColor','blue');
                    drawnow
                    hold on;

                    delete(findall(gcf,'Tag','annotation_plot'));
                    annotation('textbox',[.2 .7 .2 .2],'String',{['         Iteration   ' num2str(iii)]},...
                         'EdgeColor','white','Tag','annotation_plot','FitBoxToText','on');
                    grid off;
                    hold off;

                    subplot(8,2,[12 14 16]);
                    eqvCurvature(:) = out(:,4) - out(:,2).*gravitational + centrifugal*(out(:,1).^2);    % equivalent curvature of droplet
                    plot(1:N,eqvCurvature(:),'g', 'LineWidth',1);
                    axis([1 N -3 0]);
                    set(gca, 'XTick', [1 (N+1)/2 N], 'XTickLabel', [0 0.5 1]);
                    xlabel('\it{u}', 'Color', 'k', 'FontWeight', 'bold');
                    ylabel('\chi_{eq}          ', 'Rotation',0, 'Color', 'k', 'FontWeight', 'bold');
                    drawnow
                    hold on;
                    plot(movingInd1, eqvCurvature(movingInd1), '-s', 'MarkerEdgeColor','blue','MarkerFaceColor','blue');
                    drawnow
                    plot(movingInd2, eqvCurvature(movingInd2), '-s', 'MarkerEdgeColor','blue','MarkerFaceColor','blue');
                    drawnow

                    hold off;

                    subplot(8,2,[11 13 15]);
                    % Plot difference in volume with iteration of while loop of +PN operator
                    plot(1:iii, delVolume_arr(1:iii));
                    axis([0 iii -abs(delVolumeFirst) abs(delVolumeFirst)]);
                    drawnow;
                end
            end 
        end
        iii = iii+1;
     
    end
end



