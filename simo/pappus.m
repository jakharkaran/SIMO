function [ volume, surfaceArea, centroidVoly ]= pappus(out)
%% Calculates volume and surface area using Pappus-Guldinus theorem/Pappus Centroid theorem
% (usecase: axis symmetric figures)

% Initialize
Ax = 0;
area = 0;
Vy = 0;
N = length(out(:,1));
Px = 0;

if (rem(N,2)==0)
    % Check for even/odd number of points of point cloud 
    
    % For even number of points
    % N/2 +1 or one half + 1 of points stored in seperate arrays

    x=out(1:(N)/2+1,1);
    y=out(1:(N)/2+1,2);

   for i=1:length(x)-2
       % Calculating area under 1st to N/2 points
       Ax = Ax + ((x(i)+x(i+1))/2)*1/2*(y(i)+y(i+1))*(x(i)-x(i+1));
       area=area+1/2*(y(i)+y(i+1))*(x(i)-x(i+1));
       
       Px = Px + ((x(i)+x(i+1))/2)*sqrt(((y(i+1)-y(i))^2)+((x(i+1)-x(i))^2));
   
       Vy = Vy + (((y(i)+y(i+1))/2)/2) * (2*pi*(((y(i)+y(i+1))/2)*(x(i)-x(i+1))))*...
          (x(i)+x(i+1))/2;
   end
   
   % Calculating half of area between N/2 and N/2 + 1 point
   i = i+1;
   Ax = Ax + ((x(i)+(x(i)+x(i+1))/2)/2)*1/2*(y(i)+(y(i)+y(i+1))/2)*(x(i)-(x(i)+x(i+1))/2);
   area=area+1/2*(y(i)+(y(i)+y(i+1))/2)*(x(i)-(x(i)+x(i+1))/2)*1/2;
   
   Px = Px + ((x(i)+(x(i)+x(i+1))/2)/2) * sqrt(((y(i)-(y(i)+y(i+1))/2)^2)+((x(i)-(x(i+1)+x(i))/2)^2));
   
   Vy = Vy + (((y(i)+(y(i)+y(i+1))/2)/2)/2) * (2*pi*((y(i)+(y(i)+y(i+1))/2)/2)*(x(i)-(x(i)+x(i+1))/2)) * ...
       (x(i)+(x(i)+x(i+1))/2)/2;
else
   % For odd number of points
   % (N/2)-1 or one half of points stored in separate arrays
    
   x=out(1:((N)+1)/2,1);
   y=out(1:((N)+1)/2,2); 

   for i=1:length(x)-1
   Ax = Ax + ((x(i)+x(i+1))/2)*1/2*(y(i)+y(i+1))*(x(i)-x(i+1));  
   area = area+1/2*(y(i)+y(i+1))*(x(i)-x(i+1));
   
   Px = Px + ((x(i)+x(i+1))/2)*sqrt(((y(i+1)-y(i))^2)+((x(i+1)-x(i))^2));
   
   Vy = Vy + (((y(i)+y(i+1))/2)/2) * (2*pi*(((y(i)+y(i+1))/2)*(x(i)-x(i+1))))*...
       (x(i)+x(i+1))/2;
   end
end

% area = area of one half of point cloud
% Centroid of one half of area
centroidAreax = Ax/area;

% Perimeter of curve
perimeter = 0;
for i=1:N-1
perimeter = perimeter + sqrt((out(i+1,1)-out(i,1))^2+(out(i+1,2)-out(i,2))^2);
end

% Centroid of one half of curve
centroidCurvex = Px/(perimeter/2);

% Pappus-Guldinus theorem/Pappus Centroid theorem
surfaceArea = 2*pi*centroidCurvex*(perimeter/2);
volume = 2*pi*centroidAreax*area;

% Centroid of volume
centroidVoly = Vy/volume;
end