%% SIMO: Spline-based Interface Modeling and Optimization
% Measurement of Bo and Contact angle for a Axisymmetric drops
%% 
function [out, AlgorithmBo, Algorithmca, BoMeanArray, errorInterfaceMeanArray] = SIMOOld(pp, BoStart, BoEnd, BoIncrement)

% clear;
% close all;

ticCountDown = tic;

% plt = 1 if you want to see the plots as the code is running
plt = 0;
% captureVideo = 1 if you want to capture video of runnng plot
% plt = 1 to capture video
captureVideo = 0;

if captureVideo == 1
    filename= ['InverseAnalysis_axis_sessile', num2str(round(rand()*1000)), '.avi'];     
    vidObj = VideoWriter(filename);
    open(vidObj);
end

% Vdot*dt is the increment in volume for each iteration of +PN operator
% -PN operator will compesate this volume
Vdot = 1;
dtLarge = 0.0005;   % dt
dtSmall = 0.0001;   % Refines results with small value of dt

centrifugal = 0;    % Assuming zero cen

% Number of points on the curve
N = 50;

% Max. no. of iterations to run for each assumed value of Bo
M = 100000;

%% Digitized profile of droplet
% Input coordinates of drop profile for Bo calcualtion

% theta = 150;
% R=1;
% alpha=90-theta:2*theta/N:90+theta;  
% pp = [R*cosd(alpha') (R*sind(alpha')-R*cosd(theta))];

% csvread('file',R,C) reads data from starting at row offset R and column offset C.
% To specify the first value R,C = 0
% pp = csvread('input.csv');
% pp = edgeDetectFxn('2_1ul_10fps.png');

% BoStart = 0;      % First assumed Bo value
% BoEnd = -0.4;       % Last assumed Bo value
% BoIncrement = -0.1; % Incremental/Decremental step in assumed value of Bo as it is incremented/decremented from BoStart to BoEnd

% Calculate no. of iterations for outer loop
Bo = BoStart;
if BoIncrement == 0
    MBo = 1;
else
    MBo = ceil((BoEnd - BoStart)/BoIncrement) + 1;
end

%%
% Fitting spline over the digitized profile
outOld = mvsplint(pp,length(pp));
out = outOld;

% Normalizing drop volume, accordingly scaling drop profile coordinates 
volumeOld = pappus(out);
ReqOld = (volumeOld*3/(4*pi))^(1/3);
ReqNew = (3/(4*pi))^(1/3);
normalizingFactor = ReqNew/ReqOld;

out_first(:,1:2) = out(:,1:2) * normalizingFactor;

%%
% Smoothening drop profile
out = errorReduction(out_first);    % Makes drop profile symmetric about y axis
NN = length(out);
% Stepwise decrement of no. of coordinates in digitized drop profile to N
while NN > N
	out = mvsplint(out,NN);
	NN = NN - 1;
end
out = errorReduction(out);
%%
[outFirst,normalFirst] = mvsplint(out,N);
out = outFirst;
normalVec = normalFirst;
[volumeFirst, surfaceAreaFirst, centroidyFirst] = pappus(out);

% Characteristic dimension of droplet
Req = (volumeFirst*3/(4*pi))^(1/3);

ReqFirst = Req;
baseRadiusFirst = outFirst(1,1);
arealvFirst = surfaceAreaFirst;
areaslFirst = 2*pi*baseRadiusFirst;

%%
perimeterFirst = 0;
% Calculation of the perimeter
for i=1:N-1
perimeterFirst = perimeterFirst + sqrt((outFirst(i+1,1)-outFirst(i,1))^2+...
    (outFirst(i+1,2)-outFirst(i,2))^2);
end

% Calculate contact angles at two contact points
if out(N,3)<=0 
    lca = 180+(out(N,3));
else 
    lca=out(N,3);
end
    
if out(1,3)>=0
    rca = 180-out(1,3);
else
    rca =abs(out(1,3) );        %abs to eliminate negative - will not work for theta>90
end

caFirst = (lca+rca)/2;

%%
% Pre allocate matrices
volumeArray = zeros(MBo,N);
ReqArray = zeros(MBo,N);
capillaryConstantArray = zeros(MBo,N);
meanCurvatureArray = zeros(MBo,N,N);
eqvCurvatureArray = zeros(MBo,N,N);
stdEqvCurvatureArray = zeros(MBo,N);
meanEqvCurvatureArray = zeros(MBo,N);
coefficientOfVariationArray = zeros(MBo,N);
errorInterfaceArray = zeros(MBo,N);
perimeterArray = zeros(MBo,N);
errorPerimeterArray = zeros(MBo,N);
lcaArray = zeros(MBo,N);
rcaArray = zeros(MBo,N);
outXArray = zeros(MBo,N,N);
outYArray = zeros(MBo,N,N);
BoArray = zeros(MBo,N);
baseRadiusArray = zeros(MBo,N);
countIterationArray = zeros(MBo,1);
coefficientOfVariation = zeros(N,1);
energyArray = zeros(MBo,N);
energyArray1 = zeros(MBo,N);
energyArray2 = zeros(MBo,N);
counter = 0;

for countBo = 1:MBo
    
    % Equalizing parameters to initial value for every assumed value of Bo
	% Equivalent curvature term taking gravity into effect
    gravitational = Bo/(2*(ReqFirst^2));
    
    out = outFirst;
    normalVec = normalFirst;
    volume = volumeFirst;
    dt = dtLarge;
    
    for countIteration = 1:M
%% 
         % +PN/-PN operators
        [outNew, movingInd1, movingInd2] = pin_droplet_sync_add(out, N, volume + Vdot*dt,normalVec, gravitational, centrifugal);                       
        [out, normalVec] = mvsplint(outNew,N);

        [outNew, movingInd3, movingInd4] = pin_droplet_sync_sub(out, N, volume , normalVec, gravitational, centrifugal);
        [out, normalVec] = mvsplint(outNew,N);

        volume  = pappus(out);   
        Req = (volume*3/(4*pi))^(1/3);
        gravitational = Bo/(2*(Req^2));

        if dt == dtSmall
            % Smaller dt value for precision in measurement
            % Limit the Volume change wrt Initial volume to 0.01%
            if volumeFirst - volume >= 0.0001 * volumeFirst
                [outNew, movingInd1, movingInd2] = pin_droplet_sync_add(out, N, volume + Vdot*dt,normalVec, gravitational, centrifugal);                       
                [out, normalVec] = mvsplint(outNew,N);
            elseif volumeFirst - volume <= -0.0001 * volumeFirst
                [outNew, movingInd3, movingInd4] = pin_droplet_sync_sub(out, N, volume - Vdot*dt, normalVec, gravitational, centrifugal);
                [out, normalVec] = mvsplint(outNew,N);
            end  
        end
%%
        % Mean Curvature
        meanCurvatureArray(countBo,countIteration,:) = out(:,4);
        % Equivalent Curvature
        eqvCurvatureArray(countBo,countIteration,:) = out(:,4) - out(:,2).*gravitational;  
        eqvCurvature = out(:,4) - out(:,2).*gravitational;  
        % Standard Deviation of Equivalent Curvature
        stdEqvCurvatureArray(countBo,countIteration) = std(out(:,4) - out(:,2).*gravitational);
        % Mean of equivalent curvature
        meanEqvCurvatureArray(countBo,countIteration) = ...
            mean(eqvCurvatureArray(countBo,countIteration,:));
        % Coefficient of Variation
        coefficientOfVariationArray(countBo,countIteration) = ...
            stdEqvCurvatureArray(countBo,countIteration)./abs(meanEqvCurvatureArray(countBo,countIteration));
        coefficientOfVariation(countIteration) = abs(std(out(:,4) - out(:,2).*gravitational)/(mean(out(:,4) - out(:,2).*gravitational)));
        
        % Root Mean Square Deviation of current drop profile wrt intial drop profile
        errorInterfaceArray(countBo,countIteration) = ...
            sqrt(sum(((outFirst(:,1)-out(:,1)).^2+(outFirst(:,2)-out(:,2)).^2)/N));

        peri = 0;
        % Calculation of the perimeter
        for i=1:N-1
        peri=peri+sqrt((out(i+1,1)-out(i,1))^2+(out(i+1,2)-out(i,2))^2);
        end
        perimeterArray(countBo,countIteration) = peri;
        % Error in Perimeter of current drop profie wrt intial drop profile
        errorPerimeterArray(countBo,countIteration) = ...
            (perimeterFirst - perimeterArray(countBo,countIteration))/perimeterFirst*100;
        
        % Contact Angle
        if out(N,3)<=0
            lcaArray(countBo,countIteration)= 180+(out(N,3));
        else
            lcaArray(countBo,countIteration) = out(N,3);
        end

        if out(1,3)>=0
            rcaArray(countBo,countIteration) = 180-out(1,3);
        else
            rcaArray(countBo,countIteration) = abs(out(1,3));                        %abs to eliminate negative - will not work for theta>90
        end
        lca = lcaArray(countBo,countIteration);
        rca = rcaArray(countBo,countIteration);

        outXArray(countBo,countIteration,:) = out(:,1);
        outYArray(countBo,countIteration,:) = out(:,2);

        [volume, surfaceArea, centroidy] = pappus(out);
        Req = (volume*3/(4*pi))^(1/3);
        gravitational = Bo/(2*(Req^2));
        baseRadius = out(1,1);
        capillaryConstant = 2*gravitational;
        
        volumeArray(countBo,countIteration) = volume;
        ReqArray(countBo,countIteration) = Req;
        capillaryConstantArray(countBo,countIteration) = capillaryConstant;
        BoArray(countBo,countIteration) = Bo;
        baseRadiusArray(countBo,countIteration) = baseRadius;
        countIterationArray(countBo) = countIteration;
        
        arealv = surfaceArea;
        areasl = 2*pi*baseRadius;
        
        % Interfacial Energy
        energyArray1(countBo,countIteration) = (arealv-arealvFirst) - cosd(lca)*(areasl-areaslFirst);
        % Potential Energy
        energyArray2(countBo,countIteration) = capillaryConstant*(centroidy-centroidyFirst);
        % Total Energy
        energyArray(countBo,countIteration) = energyArray1(countBo,countIteration) + energyArray2(countBo,countIteration);

        if plt==1
            %%
            h = figure(1);
            set(h, 'Position', [50 50 1024 640], 'Color', 'white')
            
            % Plot drop shape
            subplot(8,3,[1 15]);
            plot(out_first(:,1),out_first(:,2),'k', 'LineWidth',1);
            hold on;
            plot(out(:,1),out(:,2), 'LineWidth',1);
            daspect([1,1,1]);
            axis([-2 2 0 1.5]);
            
            % Points last perturbed in +PN and -PN operator
            plot(out(movingInd1,1), out(movingInd1,2), '-s', 'MarkerEdgeColor','blue','MarkerFaceColor','blue');
            plot(out(movingInd2,1), out(movingInd2,2), '-s', 'MarkerEdgeColor','blue','MarkerFaceColor','blue');
            plot(out(movingInd3,1), out(movingInd3,2), '-o', 'MarkerEdgeColor','red','MarkerFaceColor','red');
            plot(out(movingInd4,1), out(movingInd4,2), '-o', 'MarkerEdgeColor','red','MarkerFaceColor','red');
            hold on;

            delete(findall(gcf,'Tag','annotation_plot'));

            annotation('textbox', [.15 .7 .4 .2], 'String', {['Current Bo              = ' num2str(Bo)],...
                ['Contact Angle        = ' num2str(lca)],...
                ['Error Interface        = ' num2str(errorInterfaceArray(countBo,countIteration))],...
                ['Error Perimeter      = ' num2str(errorPerimeterArray(countBo,countIteration))],...
                ['Iteration                  = ' num2str(countIteration) ],...
                ['dt = ' num2str(dt)]},...
                'EdgeColor','white','Color', 'k','Tag','annotation_plot','FitBoxToText','on');
            grid off;
            hold off;  
            
            % Plot Equivalent curvature
            subplot(8,3,[16 19 22]);
            plot(1:N,eqvCurvature,'g', 'LineWidth',1);
            set(gca, 'XTick', [1 (N+1)/2 N], 'XTickLabel', [0 0.5 1]);
            xlabel('\it{u}', 'Color', 'k', 'FontWeight', 'bold');
            ylabel('\chi_{eq}          ', 'Rotation',0, 'Color', 'k', 'FontWeight', 'bold');
            axis([1 N -2 0]);
            hold on;
            
            % Points last perturbed in +PN and -PN operator
            plot(movingInd1, eqvCurvatureArray(countBo,countIteration,movingInd1),...
                '-s', 'MarkerEdgeColor','blue','MarkerFaceColor','blue');
            plot(movingInd2, eqvCurvatureArray(countBo,countIteration,movingInd2),...
                '-s', 'MarkerEdgeColor','blue','MarkerFaceColor','blue');
            plot(movingInd3, eqvCurvatureArray(countBo,countIteration,movingInd3),...
                '-o', 'MarkerEdgeColor','red','MarkerFaceColor','red');
            plot(movingInd4, eqvCurvatureArray(countBo,countIteration,movingInd4),...
                '-o', 'MarkerEdgeColor','red','MarkerFaceColor','red');
            hold off;  
            
            % Plot Coefficient of Variation with Iteration
            subplot(8,3,[17 20 23]);
            plot(1:countIteration, coefficientOfVariationArray(countBo,1:countIteration), 'c', 'LineWidth',1);
            xlabel('Iteration #', 'Color', 'k', 'FontWeight', 'bold');
            ylabel('c_{v}     ', 'Rotation', 0, 'Color', 'k', 'FontWeight', 'bold');
            axis([0 countIteration 0 .01]); 
            
            subplot(8,3,[18 21 24]);
            plot(1:countIteration, energyArray(countBo,1:countIteration), 'c', 'LineWidth',1);
            hold on;
            plot(1:countIteration, energyArray1(countBo,1:countIteration), 'g', 'LineWidth',1);
            plot(1:countIteration, energyArray2(countBo,1:countIteration), 'r', 'LineWidth',1);
            xlabel('Iteration #', 'Color', 'k', 'FontWeight', 'bold');
            ylabel('\Delta E     ', 'Rotation', 0, 'Color', 'k', 'FontWeight', 'bold');
            axis([0 countIteration -.006 .006]); 
            hold off;
            
            % Capture frames for video
            if captureVideo == 1
               if mod(countBo,1) == 0
                   F = getframe(h);
                   writeVideo(vidObj,F);
               end
            end

        end 
            
            counter = counter + 1;
            if dt == dtLarge && countIteration>20
                % Decreasig value of dt for precision
                if mean(coefficientOfVariationArray(countBo,countIteration-9:countIteration)) >= ...
                        mean(coefficientOfVariationArray(countBo,countIteration-19:countIteration-10))
                    dt = dtSmall;
                    counter = 0;
                end
            end

            if dt == dtSmall && counter>15 
                % Breaking loop after equilibrium shape is obtained
                if mean(coefficientOfVariationArray(countBo,countIteration-4:countIteration)) >= ...
                        mean(coefficientOfVariationArray(countBo,countIteration-9:countIteration-5))
                    break;
                end
            end  
    end
Bo = Bo + BoIncrement
end

errorInterfaceMeanArray = zeros(countBo,1);
baseRadiusMeanArray = zeros(countBo,1);
volumeMeanArray = zeros(countBo,1);
lcaMeanArray = zeros(countBo,1);
errorPerimeterMeanArray = zeros(countBo,1);
BoMeanArray = zeros(countBo,1);
energyMeanArray = zeros(countBo,1);

% Taking mean of parameters for last 10 iterations for every assumed value of Bo
for countBo = 1:MBo
    
    errorInterfaceMeanArray(countBo) = mean(errorInterfaceArray(countBo,...
        countIterationArray(countBo)-9:countIterationArray(countBo)));

    baseRadiusMeanArray(countBo) = mean(baseRadiusArray(countBo,...
        countIterationArray(countBo)-9:countIterationArray(countBo)));

    volumeMeanArray(countBo) = mean(volumeArray(countBo,...
        countIterationArray(countBo)-9:countIterationArray(countBo)));

    lcaMeanArray(countBo) = mean(lcaArray(countBo,...
        countIterationArray(countBo)-9:countIterationArray(countBo)));

    errorPerimeterMeanArray(countBo) = mean(errorPerimeterArray(countBo,...
        countIterationArray(countBo)-9:countIterationArray(countBo)));

    BoMeanArray(countBo) = mean(BoArray(countBo,...
        countIterationArray(countBo)-9:countIterationArray(countBo)));
    
    energyMeanArray(countBo) = mean(energyArray(countBo,...
        countIterationArray(countBo)-9:countIterationArray(countBo)));

end

    
[errorInterrfaceMin, errorInterfaceMinInd] = min(errorInterfaceMeanArray(:));
[errorPerimeterMin, errorPerimeterMinInd] = min(abs(errorPerimeterMeanArray(:)));
   
% Bo value corresponding to minimum value of RMSD
AlgorithmBo = BoMeanArray(errorInterfaceMinInd);
Algorithmca = lcaMeanArray(errorInterfaceMinInd);
    
% Bo value corresponding to minimum error in perimeter
AlgorithmPeriBo = BoMeanArray(errorPerimeterMinInd);
AlgorithmPerica = lcaMeanArray(errorPerimeterMinInd);
    
out = [];
out(:,1) = outXArray(errorInterfaceMinInd,countIterationArray(errorInterfaceMinInd),:);
out(:,2) = outYArray(errorInterfaceMinInd,countIterationArray(errorInterfaceMinInd),:);

out = out/normalizingFactor;

if captureVideo == 1
    close(vidObj);
end

totalTime = toc(ticCountDown);