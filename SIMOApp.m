function SIMOApp

%  Create and then hide the UI as it is being constructed.
hFigure1 = figure('Name', 'SIMO ', 'MenuBar', 'none', 'ToolBar', 'none', ...
    'NumberTitle', 'off', 'Units', 'Normalized', 'OuterPosition',[0,0.05,1,0.95]);

hFig1SIMOText1 = uicontrol(hFigure1, 'Style', 'text', 'String', 'Input parameters and/or select drop profile ', ...
     'BackgroundColor', [.9 .9 .9], 'FontSize', 10, 'Units', 'Normalized', 'Position', [0 0.97 1 0.03]);
hFig1SIMOText2 = uicontrol(hFigure1, 'Style', 'text', 'String', 'SIMO', ...
    'FontSize', 70, 'Units', 'Normalized', 'Position', [0.3 0.75 0.4 0.2]);
hFig1SIMOText3 = uicontrol(hFigure1, 'Style', 'text', 'String', 'Spline-based Interface Modeling and Optimization', ...
    'FontSize', 8, 'Units', 'Normalized', 'Position', [0.3 0.7 0.4 0.1]);
hFigSIMOText4 = uicontrol(hFigure1, 'Style', 'text', 'String', 'Developed at TFTL, IIT Patna  ', ...
    'HorizontalAlignment', 'right', 'Units', 'Normalized', 'Position', [0 0 1 0.025]);

hSTPanel = uipanel(hFigure1, 'Title', 'Required for Surface Tension Calculation',...
    'Units', 'Normalized', 'Position', [0.35 0.35 0.3 0.3]); 

hdensityText1 = uicontrol(hSTPanel, 'Style', 'text', 'String', 'Density', ...
    'Units', 'Normalized', 'Position', [0.1 0.8 0.25 0.1]);
hdensityEdit = uicontrol(hSTPanel, 'Style', 'edit', 'Units', 'Normalized', ...
    'Position', [0.4 0.77 0.25 0.15], 'Callback', @hdensityEditCallback);
hdensityText2 = uicontrol(hSTPanel, 'Style', 'text', 'String', 'kg/m3', ...
    'Units', 'Normalized', 'Position', [0.7 0.8 0.25 0.1]);

hSTButtonGroup = uibuttongroup(hSTPanel, 'Units','Normalized', 'Title', 'Input one of the Following',...
    'Position',[0 0 1 0.7], 'Selection', 'on', 'SelectionChangedFcn',{@hSTButtonGroupCallback});
              
% Create three radio buttons in the button group.
hVolumeRadioButton = uicontrol(hSTButtonGroup,'Style', 'radiobutton','String','Volume',...
                  'Units','Normalized', 'Position',[0.1 0.7 0.25 0.2]);

hVolumeEdit = uicontrol(hSTButtonGroup,'Style', 'edit', 'Units','Normalized', ...
    'Position',[0.4 0.65 0.25 0.25], 'Callback', @hVolumeEditCallback);
hVolumePopUpMenu = uicontrol(hSTButtonGroup,'Style', 'popupmenu', ...
    'String', {'cubic metre', 'microlitre', 'nanolitre'}, ...
    'Units','Normalized', 'Position',[0.7 0.65 0.25 0.25]);
           
hScaleRadioButton = uicontrol(hSTButtonGroup,'Style', 'radiobutton','String','Scale',...
                  'Units','Normalized', 'Position',[0.1 0.4 0.25 0.2]);
hScaleEdit = uicontrol(hSTButtonGroup,'Style', 'edit', 'Units','Normalized', 'Visible', 'off',...
    'Position',[0.4 0.35 0.25 0.25], 'Callback', @hScaleEditCallback);
hScalePopUpMenu = uicontrol(hSTButtonGroup,'Style', 'popupmenu', 'Visible', 'off', 'String', ...
    {'pixel per milimetre', 'pixel per metre'}, 'Units','Normalized', 'Position',[0.7 0.35 0.25 0.25],....
    'Callback', @hScalePopUpMenuCallback);
hScaleInImage = uicontrol(hSTButtonGroup, 'Style', 'pushbutton', 'String', 'Select Scale from Image',...
    'Visible', 'off', 'Units','Normalized', 'Position', [0.4 0.05 0.55 0.25], ...
    'Callback', @hScaleInImageCallback);

hScalePanel = uipanel(hSTPanel, 'Units','Normalized', 'Visible', 'off',...
    'Title', 'Input Scale Conversion', 'Position',[0 0 1 0.7]);
hScaleInImageText = uicontrol(hScalePanel, 'Style', 'text', 'String', 'Scale', ...
    'Units','Normalized', 'Position',[0.1 0.35 0.25 0.25]);
hScaleInImageEdit = uicontrol(hScalePanel,'Style', 'edit', 'Units','Normalized',...
    'Position',[0.4 0.4 0.25 0.25], 'Callback', @hScaleInImageEditCallback);
hScaleInImagePopUpMenu = uicontrol(hScalePanel,'Style', 'popupmenu', 'String', ...
    {'milimetre', 'metre'}, 'Units','Normalized', 'Position',[0.7 0.4 0.25 0.25],....
    'Callback', @hScaleInImagePopUpMenuCallback);

hSelectProfile = uicontrol(hFigure1, 'Style', 'pushbutton', 'String', 'Select Drop Profile',...
    'Units','Normalized', 'Position', [0.35 0.15 0.3 0.075], ...
    'Callback', @hSelectProfileCallback);

%%

rectImageExist = 0;
finalPosition = [];
ImageOrig = [];
ImageRotated = [];
ImageCrop = [];
ImageCropCanny = [];
pp = [];
boundary = [];
scaleLength = [];


%%
hFigure2 = figure('Name', 'SIMO ', 'Visible', 'off', 'MenuBar', 'none', ...
    'ToolBar', 'none', 'NumberTitle', 'off', 'Units', 'Normalized', 'OuterPosition',[0,0.05,1,0.95]);
            
% hFigure1.Visible = 'on';
% hFigure2.Visible = 'off';

% ImageOrig = imread('1.bmp');
% ImageOrigSize = size(ImageOrig);
% ImageOrigDimension = [num2str(ImageOrigSize(2)) ' X ' num2str(ImageOrigSize(1)) ' pixels'];

hFig2SIMOText1 = uicontrol(hFigure2, 'Style', 'text', 'String', 'Rotate Image for Contact Surface to face upwards', ...
    'BackgroundColor', [.9 .9 .9], 'FontSize', 10, 'Units', 'Normalized', 'Position', [0 0.97 1 0.03]);
hImageRotatedDimText = uicontrol(hFigure2, 'Style', 'text', 'HorizontalAlignment', 'left', ...
    'Units', 'Normalized', 'Position', [0.01 0.925 0.09 0.025]);
hImageCropDimText = uicontrol(hFigure2, 'Style', 'text', 'HorizontalAlignment', 'left',...
    'Visible', 'off', 'Units', 'Normalized', 'Position', [0.505 0.925 0.11 0.025]);
hVolumeUnScaledText = uicontrol(hFigure2, 'Style', 'text', 'HorizontalAlignment', 'left',...
    'Visible', 'off', 'Units', 'Normalized', 'Position', [0.505 0.9 0.11 0.025]);
hVolumeScaledText = uicontrol(hFigure2, 'Style', 'text', 'HorizontalAlignment', 'left',...
    'Visible', 'off', 'Units', 'Normalized', 'Position', [0.505 0.876 0.11 0.025]);

% Construct the components.
hImageOrig = axes('Parent', hFigure2, 'Units','Normalized','Position',[0.01,0.3,0.485,0.65]); 
hImageCrop = axes('Parent', hFigure2, 'Units','Normalized','Position',[0.505,0.3,0.485,0.65]);
axis off;
% hEdgeDetectCanny = axes('Parent', hFigure2, 'Units','Pixels','Position',[660,10,173,130]);
% axis off;
hEdgeDetectProfile = axes('Parent', hFigure2, 'Units','Normalized', 'Position',[0.84,0.8,0.15,0.15]);
axis off;

hRotatePanel = uipanel(hFigure2, 'BorderType', 'None', 'Units', 'Normalized',...
    'Position', [0.01 0.25 0.485 0.045]); 
hRotateEdit = uicontrol(hRotatePanel, 'Style', 'edit', 'Units', 'Normalized', 'String', '0',...
    'Position', [0.6 0 .15 1], 'Callback', {@hRotateEditCallback});
hRotate180Button = uicontrol(hRotatePanel, 'Style', 'pushbutton', 'Units', 'Normalized', ...
    'String', 'Rotate by 180 degrees', 'Position', [0.1 0 .2 1],...
    'Callback', {@hRotate180ButtonCallback, hRotateEdit});
hRotateCounterClockwiseButton = uicontrol(hRotatePanel, 'Style', 'pushbutton', 'Units', 'Normalized', ...
    'String', 'Rotate Counter Clockwise', 'Position', [0.4 0 .15 1], ...
    'Callback', {@hRotateCounterClockwiseButtonCallback, hRotateEdit});
hRotateClockwiseButton = uicontrol(hRotatePanel, 'Style', 'pushbutton', 'Units', 'Normalized', ...
    'String', 'Rotate Clockwise', 'Position', [0.8 0 .15 1], ...
    'Callback', {@hRotateClockwiseButtonCallback, hRotateEdit});

hrectImage = uicontrol(hFigure2, 'Style', 'pushbutton', 'Units', 'Normalized',...
    'String', 'Contact Surface is Facing Upwards', ...
    'Position', [0.03 0.15 0.445 0.05], 'Callback', @hrectImageCallback);

hCompute = uicontrol(hFigure2, 'Style', 'pushbutton', 'Units', 'Normalized',...
    'String', 'I am satisfied! Compute.', 'Visible', 'off',...
    'Position', [0.555,0.05,0.385,0.07], 'Callback', @hComputeCallback);

hImageCropThresholdPanel = uipanel(hFigure2, 'Visible', 'off', ...
    'BorderType', 'None', 'Units', 'Normalized', 'Title', 'Threshold', ...
    'TitlePosition', 'centertop', 'Position', [0.05 0.05 0.15 0.08]); 
hImageCropThreshold = uicontrol(hImageCropThresholdPanel, 'Style','slider', 'min', 0, 'max', 0.9999999,...
    'Value', 0.5, 'Units', 'Normalized', 'Position',[0.1,0.1,0.8,0.8], ...
    'Callback', @hImageCropThresholdCallback);
hImageCropThresholdTextMin = uicontrol(hImageCropThresholdPanel, 'Style','text', ...
    'String','1', 'Units', 'Normalized', 'Position',[0.93,0.1,0.05,0.6]);
hImageCropThresholdTextMax = uicontrol(hImageCropThresholdPanel, 'Style','text', ...
    'String','0', 'Units', 'Normalized', 'Position',[0.02,0.1,0.05,0.6]);

hBoRangePanel = uipanel(hFigure2, 'Visible', 'off', 'Units', 'Normalized', 'Title', 'Bond number Range',...
    'Position', [0.555, 0.15, 0.385, 0.07]);
hBoStartText = uicontrol(hBoRangePanel, 'Style', 'text', 'Units', 'Normalized',...
    'String', 'Bo Start', 'Position', [0 0 0.25 0.8]);
hBoStartEdit = uicontrol(hBoRangePanel, 'Style', 'edit', 'Units', 'Normalized',...
    'String', '-2', 'Position', [0.2 0.1 0.25 0.9], 'Callback', @hBoStartEditCallback);
hBoEndText = uicontrol(hBoRangePanel, 'Style', 'text', 'Units', 'Normalized',...
    'String', 'Bo End', 'Position', [0.5 0 0.25 0.8]);
hBoEndEdit = uicontrol(hBoRangePanel, 'Style', 'edit', 'Units', 'Normalized',...
    'String', '10', 'Position', [0.7 0.1 0.25 0.9], 'Callback', @hBoEndEditCallback);

%%

BoStart = str2num(hBoStartEdit.String);
BoEnd = str2num(hBoEndEdit.String);
% ImageOrig = imread([pathname,filename]);
% imshow(ImageOrig, 'Parent', hImageOrig);
% imshow(ImageCrop, 'Parent', hImageCrop);
% imshow([pathname,filename], 'Parent', hImageOrig)

%%
    function hdensityEditCallback(source,~)
        if isempty(source.String)
            
        elseif isempty(str2num(source.String))
            source.String = [];
            warndlg('Density input must be numerical');
        end  
        
    end

    function hSTButtonGroupCallback(~, event)
        if strcmp(event.NewValue.String,'Volume')
            hScaleEdit.String = [];
            hScaleEdit.Visible = 'off';
            hScalePopUpMenu.Visible = 'off';
            hVolumeEdit.Visible = 'on';
            hVolumePopUpMenu.Visible = 'on'; 
            hVolumeScaledText.Visible= 'off';
            if strcmp(hFigure2.Visible, 'on')
                hScaleInImage.Visible = 'off';
            end
        elseif strcmp(event.NewValue.String,'Scale')
            hVolumeEdit.String = [];
            hScaleEdit.Visible = 'on';
            hScalePopUpMenu.Visible = 'on';
        	hVolumeEdit.Visible = 'off';
            hVolumePopUpMenu.Visible = 'off'; 
            if strcmp(hFigure2.Visible, 'on')
                hScaleInImage.Visible = 'on';
            end
        end
    end

    function hVolumeEditCallback(source, ~)
        if isempty(source.String)
            
        elseif isempty(str2num(source.String))
            source.String = [];
            warndlg('Volume input must be numerical');
        end   
    end

    function hScaleEditCallback(source, ~)
        if isempty(source.String)
            
        elseif isempty(str2num(source.String))
            source.String = [];
            warndlg('Scale input must be numerical');
        end
        
        if (isempty(str2num(hScaleEdit.String)) ~= 1) && (isempty(pp) ~= 1)
            edgeDetect(hImageCropThreshold, 'drop');
        end 
    end 

    function hScalePopUpMenuCallback(~,~)
        if (isempty(str2num(hScaleEdit.String)) ~= 1) && (isempty(pp) ~= 1)
            edgeDetect(hImageCropThreshold, 'drop'); 
        end
    end

    function hSelectProfileCallback(~, ~)

        [filename, pathname] = uigetfile({'*.jpg'; '*.bmp'; '*.png'; '*.tiff'}, 'Select a Drop Profile');
        if ischar(filename) && ischar(pathname)
            hFigure1.Visible = 'off';
            hFigure2.Visible = 'on';
            
            filepath = strcat(pathname,filename);
            image = imread(filepath);
            
            if size(image,3) == 4
                ImageOrig = image(:,:,1:3);
            else
                ImageOrig = image;
            end
            
            ImageRotated = ImageOrig;
            ImageOrigSize = size(ImageOrig);
            hImageRotatedDimText.String = [num2str(ImageOrigSize(2)) ' X ' num2str(ImageOrigSize(1)) ' pixels'];
            imshow(ImageOrig, 'Parent', hImageOrig);
            
            if strcmp(hScalePopUpMenu.Visible, 'on')
                hScaleInImage.Visible = 'on';
            end
        end      
    end

%%
    function hScaleInImageEditCallback(source,~)
        if isempty(source.String)
            warndlg('Scale conversion requires input');
        elseif isempty(str2num(source.String))
            source.String = [];
            warndlg('Scale conversion input must be numerical');
        end 
        if isempty(str2num(hScaleInImageEdit.String)) ~= 1
            edgeDetect(hImageCropThreshold, 'scale'); 
        end
    end

    function hScaleInImagePopUpMenuCallback
        if isempty(str2num(hScaleInImageEdit.String)) ~= 1
            edgeDetect(hImageCropThreshold, 'scale'); 
        end
    end

    function hScaleInImageCallback(~,~)
        hFig2SIMOText1.String = 'Select Scale Region';
        hSTButtonGroup.Visible = 'off';
        hVolumeEdit.String = [];
        hScaleEdit.String = [];
        
        rectImageCrop = imrect(hImageOrig);
        initialPosition = getPosition(rectImageCrop);
        rectImage(rectImageCrop, initialPosition, ImageRotated, hImageCropThreshold, 'scale');
        
        hFig2SIMOText1.String = 'Move, Resize Scale Region and Input scale conversion.';

    end

%%
    function hRotate180ButtonCallback(~,~, handleHRotateEdit)
        cla(hImageCrop); 
        cla(hEdgeDetectProfile); 
        hEdgeDetectProfile.Visible = 'off';
        
        initialAngle = str2num(handleHRotateEdit.String);
        angle = mod(initialAngle + 180, 360);
        handleHRotateEdit.String = angle;
        ImageRotated = imrotate(ImageOrig, angle, 'bicubic', 'loose');
        
        imshow(ImageRotated, 'Parent', hImageOrig);
        
        if rectImageExist == 1
            rectImageCrop = imrect(hImageOrig, finalPosition);
            rectImage(rectImageCrop, finalPosition, ImageRotated, hImageCropThreshold, 'drop');
        end
    end
 
    function hRotateClockwiseButtonCallback(~,~, handleHRotateEdit)
        cla(hImageCrop); 
        cla(hEdgeDetectProfile); 
        hEdgeDetectProfile.Visible = 'off';

        initialAngle = str2num(handleHRotateEdit.String);
        angle = mod(initialAngle - 0.1, 360);
        handleHRotateEdit.String = angle;
        ImageRotated = imrotate(ImageOrig, angle, 'bicubic', 'loose');  
        ImageRotatedSize = size(ImageRotated);
        hImageRotatedDimText.String = [num2str(ImageRotatedSize(2)) ' X ' num2str(ImageRotatedSize(1)) ' pixels'];
        
        imshow(ImageRotated, 'Parent', hImageOrig);
        
        if rectImageExist == 1
            rectImageCrop = imrect(hImageOrig, finalPosition);
            rectImage(rectImageCrop, finalPosition, ImageRotated, hImageCropThreshold, 'drop');
        end
            
     end
 
    function hRotateCounterClockwiseButtonCallback(~,~, handleHRotateEdit)
        cla(hImageCrop); 
        cla(hEdgeDetectProfile); 
        hEdgeDetectProfile.Visible = 'off';
        
        initialAngle = str2num(handleHRotateEdit.String);
        angle = mod(initialAngle + 0.1, 360);
        handleHRotateEdit.String = angle;
        ImageRotated = imrotate(ImageOrig, angle, 'bicubic', 'loose');
        ImageRotatedSize = size(ImageRotated);
        hImageRotatedDimText.String = [num2str(ImageRotatedSize(2)) ' X ' num2str(ImageRotatedSize(1)) ' pixels'];
        
        imshow(ImageRotated, 'Parent', hImageOrig);
        
        if rectImageExist == 1
            rectImageCrop = imrect(hImageOrig, finalPosition);
            rectImage(rectImageCrop, finalPosition, ImageRotated, hImageCropThreshold, 'drop');
        end
    end
 
    function hRotateEditCallback(source,~) 
        if isempty(source.String)
            warndlg('Angle of rotation requires numerical input');
        elseif isempty(str2num(source.String))
            source.String = [];
            warndlg('Angle of rotation input must be numerical');
        end 
        
        cla(hImageCrop); 
        cla(hEdgeDetectProfile); 
        hEdgeDetectProfile.Visible = 'off';
        
        angle = str2num(source.String);
        ImageRotated = imrotate(ImageOrig, angle, 'bicubic', 'loose');
        ImageRotatedSize = size(ImageRotated);
        hImageRotatedDimText.String = [num2str(ImageRotatedSize(2)) ' X ' num2str(ImageRotatedSize(1)) ' pixels'];
  
        imshow(ImageRotated, 'Parent', hImageOrig);
        
        if rectImageExist == 1
            rectImageCrop = imrect(hImageOrig, finalPosition);
            rectImage(rectImageCrop, finalPosition, ImageRotated, hImageCropThreshold, 'drop');
        end
    end

%%
    function hrectImageCallback(~, ~) 
        hFig2SIMOText1.String = 'Select Drop Region';
        hFig2SIMOText1.Visible = 'on';
        hrectImage.Visible = 'off';
        
        rectImageCrop = imrect(hImageOrig);
        initialPosition = getPosition(rectImageCrop);
        rectImage(rectImageCrop, initialPosition, ImageRotated, hImageCropThreshold, 'drop');
       
        hFig2SIMOText1.String = 'Move, Resize Drop Region and Rotate Image to your satisfaction';
        hImageCropThresholdPanel.Visible = 'on';
        hSTPanel.Position = [0.25 0.01 0.25 0.23];
        hSTPanel.Parent = hFigure2;
        hCompute.Visible = 'on';
        hBoRangePanel.Visible = 'on';
        rectImageExist = 1;
        
    end

    function rectImage(rectImageCrop, pos, handleImageOrig, handlehImageCropThreshold, rectType)
           
        rectImageCrop.addNewPositionCallback(@(pos)rectImageCropCallback(pos, handleImageOrig, handlehImageCropThreshold, rectType));
        rectImageCropCallback(pos, handleImageOrig, handlehImageCropThreshold, rectType);
        
        finalPosition = getPosition(rectImageCrop);
    end

    function rectImageCropCallback(position, handleImageOrig, handlehImageCropThreshold, rectType)
        x1 = position(1);
        y1 = position(2);
        x2 = position(1) + position(3);
        y2 = position(2) + position(4);
        
        ImageCrop = handleImageOrig(round(y1:y2), round(x1:x2), :);
       
        edgeDetect(handlehImageCropThreshold, rectType);
    end

    function hImageCropThresholdCallback(source, ~)
        edgeDetect(source, 'drop');
    end

%%
    function edgeDetect(handlehImageCropThreshold, rectType)
        
        if size(ImageCrop,3) == 3
            ImageCropGray = rgb2gray(ImageCrop);
        else
            ImageCropGray = ImageCrop;
        end
        BwImageCrop = edge(ImageCropGray,'Canny', handlehImageCropThreshold.Value);
       
        switch rectType
            case 'drop'
                clearvars boundary pp;

                row = length(BwImageCrop(:,1)) - 1;

        %         imshow(BwImageCrop, 'Parent', hEdgeDetectCanny);

                [~,col] = find(BwImageCrop(row,:));
                colStart = min(col);
                colEnd = max(col);
                BwImageCrop(row,colStart:colEnd) = 1;

                boundary = bwtraceboundary(BwImageCrop,[row, colStart],'W');
                pp(:,1) = boundary(find(boundary(:,1)<row),2); %#ok<*FNDSB>
                pp(:,2) = boundary(find(boundary(:,1)<row),1);
                pp(:,2) = max(pp(:,2)) - pp(:,2);
                pp = flipud(pp);

                pp = errorReduction(pp);
                plot(pp(:,1), pp(:,2), '.', 'Parent', hEdgeDetectProfile);
                daspect(hEdgeDetectProfile, [1 1 1]);

                hVolumeUnScaledText.String = [num2str(round(pappus(pp),2)) ' cubic pixels'];
                hVolumeUnScaledText.Visible = 'on';

                ImageCropSize = size(ImageCrop);
                hImageCropDimText.String = [num2str(ImageCropSize(2)) ' X ' num2str(ImageCropSize(1)) ' pixels'];
                hImageCropDimText.Visible = 'on';

                if isempty(str2num(hScaleEdit.String)) ~= 1
                    scale = str2num(hScaleEdit.String);
                    ppScaled = pp/scale;
                    volumeScaled = pappus(ppScaled);

                    switch hScalePopUpMenu.Value
                        case 1
                            hVolumeScaledText.String = [num2str(round(volumeScaled,2)) ' cubic milimetre'];
                        case 2  
                            hVolumeScaledText.String = [num2str(round(volumeScaled,2)) ' cubic metre'];
                    end
                    hVolumeScaledText.Visible = 'on';    
                end            

                cla(hImageCrop);
                imshow(ImageCrop, 'Parent', hImageCrop); 
                hold(hImageCrop, 'on');
                plot(hImageCrop, boundary(:,2), boundary(:,1), '.r');
                hold(hImageCrop, 'off');
        
            case 'scale'
                hVolumeScaledText.Visible = 'off';
                hScalePanel.Visible = 'on';
                                
                [~,col] = find(BwImageCrop);                
                scaleLength = abs(max(col) - min(col));  
                hVolumeUnScaledText.String = [num2str(scaleLength) ' pixels'];
                
                cla(hImageCrop);
                imshow(ImageCrop, 'Parent', hImageCrop); 
                       
        end
    end

    function hComputeCallback(~,~)
        hFigure2.Visible = 'off';
        drawnow; 
        csvwrite('input.csv', pp);
        
        if isempty(str2num(hVolumeEdit.String)) ~= 1
            switch hVolumePopUpMenu.Value
                case 1
                    volumeScale = 'cubic metre';
                case 2  
                    volumeScale = 'microlitre';
                case 3
                    volumeScale = 'nanolitre';
            end
            SIMOAppCompute(pp, str2num(hdensityEdit.String), str2num(hVolumeEdit.String), volumeScale, 0, num2str(hVolumeUnScaledText.String), BoStart, BoEnd)
        
        elseif isempty(str2num(hScaleEdit.String)) ~= 1
            switch hScalePopUpMenu.Value
                case 1
                    scale = 'pixel per milimetre';
                case 2
                    scale = 'pixel per metre';
            end
            SIMOAppCompute(pp, str2num(hdensityEdit.String), str2num(hScaleEdit.String), scale, 0, num2str(hVolumeUnScaledText.String), BoStart, BoEnd)

        elseif isempty(str2num(hScaleInImageEdit.String)) ~= 1
            switch hScaleInImagePopUpMenu.Value
                case 1.
                    scale = 'milimetre';
                case 2
                    scale = 'metre';
            end
            SIMOAppCompute(pp, str2num(hdensityEdit.String), scaleLength, scale, str2num(hScaleInImageEdit.String), num2str(hVolumeUnScaledText.String), BoStart, BoEnd)

        else 
            SIMOAppCompute(pp, str2num(hdensityEdit.String), scaleLength, 0, str2num(hScaleInImageEdit.String), num2str(hVolumeUnScaledText.String), BoStart, BoEnd)

        end  
        
    end

    function hBoStartEditCallback(source,~)
        if isempty(source.String)
            
        elseif isempty(str2num(source.String))
            source.String = [];
            warndlg('Bo Start input must be numerical');
        end         
        BoStart = str2num(source.String);
    end

    function hBoEndEditCallback(source,~)
        if isempty(source.String)
            
        elseif isempty(str2num(source.String))
            source.String = [];
            warndlg('Bo End input must be numerical');
        end         
        BoEnd = str2num(source.String);
    end

 end