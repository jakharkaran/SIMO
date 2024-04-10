function [outNew, movingInd3, movingInd4, minIndex] = pin_droplet_sync_sub(outNew, N, volumeOrig, normalVec, gravitational, centrifugal)
%%  -PN operator
%    Subtracts volume by moving points of minimum curvature in the
%    direction of normal vector, decresing the volume
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

    movingInd3 = 1;
    movingInd4 = 1;
    minInd = 1;

    % Initialise no. of loops    
    iii = 1;
    %% Compensate volume addition 
    while abs(delVolume) > errTol 

        % Index of point with minimum curvature
        [~, minInd] = min(out(2:N-1,4) - out(2:N-1,2)*gravitational + centrifugal*(out(2:N-1,1).^2));
        % End points are constrained
        minInd = minInd + 1;
        
        % Calculating appropriate direction of normal vector of point with
        % minimum curvature for perturbation
        [ normalVec ] = normalVector( out, normalVec, N, volumeCurr, delVolumeFirst, minInd );
        normalVecUpdated = normalVec;

        % For axisymmetric curve two points symmetric to axis of symmetry are perturbed
        % Perturbation vector
        perturbVec = normalVec(minInd, 1:2); 
        perturbVecSync = normalVec(N+1-minInd, 1:2); 
        % Normalise perturbation vector  
        perturbVec = perturbVec/norm(perturbVec);  
        perturbVecSync = perturbVecSync/norm(perturbVecSync);

        % Perturb points with minimum curvature to decrease volume
        out(minInd,1:2) = out(minInd,1:2) - (P*(delVolume/2)*perturbVec ...
            + I*(errInt/2)*perturbVec...
            + D*(delVolume - errPrev)/2*perturbVec);

        out(N-minInd+1,1:2) = out(N-minInd+1,1:2) - (P*(delVolume)/2*perturbVecSync...
            + I*(errInt)/2*perturbVecSync...
            + D*(delVolume - errPrev)/2*perturbVecSync);

        % Storing previous volume difference
        errPrev = delVolume;
        
        % Update points with a new spline
        [outUpdated,normalVec] = mvsplint(out,N);    
        out = outUpdated;
        % Calculate the current volume
        volumeCurr=pappus(outUpdated);
        
        % Calculate the new volume difference
        delVolume = volumeOrig - volumeCurr;
        delVolumeArr(iii) = delVolume;
        
        % Compute the error integral term
        errInt = errInt + (delVolume);
        % Initialize the integral term if the error is zero
        if abs(errInt)>abs(2*volumeCurr)
            errInt = 0;
        end
        minIndex(iii) = minInd;
        %% Output
        
        outNew = out;

        % Index of perturbing point
        movingInd3 = minInd;
        movingInd4 = N-minInd+1; 
        
        %% Plot    

        if plt == 1
            if iii > 10
                if mod(iii,20)==0    
                    h = figure(1);
                    subplot(8,2,[1 10]);
                    set(h, 'Position', [50 50 1024 640], 'Color', 'white')
                    plot(out(:,1), out(:,2), 'LineWidth',1.5);
                    axis([-2.5 2.5 0 2]);
                    daspect([1,1,1]);
                    drawnow
                    hold on;
                    
                    quiver(out(minInd,1),out(minInd,2),0.1*normalVecUpdated(minInd,1), .1*normalVecUpdated(minInd,2));
                    quiver(out(N+1-minInd,1),out(N+1-minInd,2),0.1*normalVecUpdated(N+1-minInd,1), .1*normalVecUpdated(N+1-minInd,2));

                    movingInd3 = minInd;
                    movingInd4 = N-minInd+1;            

                    plot(out(movingInd3,1), out(movingInd3,2), '-o', 'MarkerEdgeColor','red','MarkerFaceColor','red');
                    drawnow
                    plot(out(movingInd4,1), out(movingInd4,2), '-o', 'MarkerEdgeColor','red','MarkerFaceColor','red');
                    drawnow
                    hold on;

                    delete(findall(gcf,'Tag','annotation_plot'));
                    annotation('textbox',[.2 .7 .2 .2],'String',{['         Iteration   ' num2str(iii)]},...
                         'EdgeColor','white','Tag','annotation_plot','FitBoxToText','on');
                    grid off;
                    hold off;
                    delVolume;
                    subplot(8,2,[12 14 16]);
                    eqvCurvature(:) = out(:,4) - out(:,2).*gravitational;    % equivalent curvature of droplet
                    plot(1:N,eqvCurvature(:),'g', 'LineWidth',1);
                    axis([1 N -3 0]);
                    set(gca, 'XTick', [1 (N+1)/2 N], 'XTickLabel', [0 0.5 1]);
                    xlabel('\it{u}', 'Color', 'k', 'FontWeight', 'bold');
                    ylabel('\chi_{eq}          ', 'Rotation',0, 'Color', 'k', 'FontWeight', 'bold');
                    drawnow
                    hold on;
                    plot(movingInd4, eqvCurvature(movingInd3), '-o', 'MarkerEdgeColor','red','MarkerFaceColor','red');
                    drawnow
                    plot(movingInd3, eqvCurvature(movingInd4), '-o', 'MarkerEdgeColor','red','MarkerFaceColor','red');
                    drawnow
                    hold off;

                    subplot(8,2,[11 13 15]);
                    % Plot difference in volume with iteration of while loop of +PN operator
                    plot(1:iii, delVolumeArr(1:iii));
                    axis([0 iii -abs(delVolumeFirst) abs(delVolumeFirst)]);
                    drawnow;
                end
            end
        end 
        iii = iii+1;
     
    end
end