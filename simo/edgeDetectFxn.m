function [ out ] = edgeDetectFxn( image )
%% Edge Detection Function
% Digitze dropet image using canny edge detector

[I2] = imread(image);
I = rgb2gray(I2);

Bw = edge(I,'Canny', .5);

% imshow(I);
% hold on;
row = 90;

[~,col] = find(Bw(row,:));
colStart = min(col);
colEnd = max(col);

Bw(1:row,:) = 0;
Bw(row,colStart:colEnd) = 1;

% imshow(Bw);
% hold on;

boundary = bwtraceboundary(Bw,[row, colStart],'E');

pp(:,1) = boundary(find(boundary(:,1)>row),2);
pp(:,2) = boundary(find(boundary(:,1)>row),1);

% plot(boundary(:,2),boundary(:,1));
% hold on;
% plot(boundary(1,2),boundary(1,1), '*');

pp = errorReduction(pp);
% plot(pp(:,1), pp(:,2), '.');
% daspect([1 1 1]);
% hold on;
% plot(pp(1,1), pp(1,2), '*');

out = pp;

end

