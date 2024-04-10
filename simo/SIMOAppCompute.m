function SIMOAppCompute( pp, density, scalingQuantity , scale, scalingQuantityPer, volumeUnScaled, BoStart, BoEnd )
%SIMOAppCompute Computes and displays results of SIMO application
% 
hFigure3 = figure('Name', 'SIMO ', 'MenuBar', 'none', 'ToolBar', 'none', ...
    'NumberTitle', 'off', 'Units', 'Normalized', 'OuterPosition',[0,0.05,1,0.95]);

hFig3SIMOText1 = uicontrol(hFigure3, 'Style', 'text', 'String', 'Computing... ', ...
     'FontSize', 50, 'Units', 'Normalized', 'Position', [0 0.5 1 0.12]);
 
drawnow;
 
[outUnScaled, AlgorithmBo, Algorithmca, BoMeanArray, errorInterfaceMeanArray] = SIMO(pp, BoStart, BoEnd);

% BoIncrement = 0.1;
% [outUnScaled, AlgorithmBo, Algorithmca, BoMeanArray, errorInterfaceMeanArray] = SIMOOld(pp, BoStart, BoEnd, BoIncrement);

SIMOdata{1,1} = 'Bond Number';
SIMOdata{2,1} = AlgorithmBo;

SIMOdata{1,2} = 'Contact Angle';
SIMOdata{2,2} = Algorithmca;

if isempty(density) ~= 1
    [outScaled, surfaceTension, volume, Req  ] = surfaceTensionFind( pp, AlgorithmBo, density, scalingQuantity, scale, scalingQuantityPer );
    
    SIMOdata{1,4} = 'Surface Tension (Newton per metre)';
    SIMOdata{2,4} = surfaceTension;
    
    SIMOdata{1,5} = 'Volume (cubic metre)';
    SIMOdata{2,5} = volume;
    
    SIMOdata{1,7} = 'Characteristic Dimension of Drop (metre)';
    SIMOdata{2,7} = Req;
    
    SIMOdata{1,9} = 'Drop Profile X coordinate (metres)';
    SIMOdata{1,10} = 'Drop Profile Y corrdinate (metres)';
    SIMOdata{2, 9} = outScaled(:,1);
    SIMOdata{2, 10} = outScaled(:,2);

end  

SIMOdata{1,6} = 'Volume (cubic pixels)';
SIMOdata{2,6} = volumeUnScaled;

SIMOdata{1,11} = 'Drop Profile X coordinate (pixels)';
SIMOdata{1,12} = 'Drop Profile Y corrdinate (pixels)';
SIMOdata{2, 11} = outUnScaled(:,1);
SIMOdata{2, 12} = outUnScaled(:,2);

close(hFigure3);
clearvars -except SIMOdata

save('SIMO.mat');


end

