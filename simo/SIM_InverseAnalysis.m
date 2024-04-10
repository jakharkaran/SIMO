%% SIM: Spline-based Interface Modeling - Inverse Analysis
% Determining Bo or Bw or both for a given digitized droplet shape by Inverse Analysis
close all;
clear;

%% Digitized profile of droplet
% Input coordinates of drop profile for Bo calcualtion

% theta = 90;
% R=1;
% alpha=90-theta:2*theta/N:90+theta;  
% pp = [R*cosd(alpha') (R*sind(alpha')-R*cosd(theta))];

% csvread('file',R,C) reads data from starting at row offset R and column offset C.
% To specify the first value R,C = 0
% pp = csvread('Bo-0.3_120_SS_1max%4.csv',0,0);
pp = edgeDetectFxn('2_1ul_10fps.png');

changeBo = 1;   % changeBo = 1 to guess values of Bo
changeBw  = 0;  % changeBw = 1 to guess values of Bw
changeBoBw = 0; % changeBoBw = 1 to guess values of both Bo and Bw

% Number of points on the curve
N = 50;

if changeBo == 1 || changeBoBw == 1
    % Range of Bo
    BoStart = 5;        % First assumed Bo value
    BoEnd = -5;         % Last assumed Bo value
    BoIncrement = -.1; % Incremental/Decremental step in assumed value of Bo as it is incremented/decremented from BoStart to BoEnd 
    Bo = BoStart;
    
    % Define value of Bw if changeBo = 1
    if changeBo == 1
        Bw = 0;
    end
end

if changeBw == 1 || changeBoBw == 1
    % Range of Bw
    BwStart = -2;       % First assumed Bw value
    BwEnd = 2;          % Last assumed Bw value
    BwIncrement = .1;   % Incremental/Decremental step in assumed value of Bw as it is incremented/decremented from BwStart to BwEnd 
    Bw = BwStart;
 
    % Define value of Bo if changeBw = 1
    if changeBw == 1
        Bo = 0;
    end
end

% Calculating No. of iterations to run for the following loop
if changeBo == 1 || changeBoBw == 1
    M1 = ceil((BoEnd - BoStart)/BoIncrement) +1;
end
if changeBw ==1 || changeBoBw == 1
    M2 = ceil((BwEnd - BwStart)/BwIncrement) +1;
end

%%
% % Fitting spline over the digitised profile
% [out,normal_vec]=mvsplint(pp,N);
% 
% % % Calculating original volume using pappus centroid theorem
% volumeOld1 = pappus(out);
% ReqOld1 = (volumeOld1*3/(4*pi))^(1/3);
% ReqNew1 = (3/(4*pi))^(1/3);
% noramalisingFactor1 = ReqNew1/ReqOld1;
% 
% out_first1(:,1:2) = out(:,1:2) * noramalisingFactor1;
% [out,~] = mvsplint(out_first1,N);
% volumeOld = pappus(out);
% % 
% % % Introducing systematic error in theoritical droplet shape    
% % 
% perturbVec = -1 + 2*rand([N-2,2]);
% errorFactor = volumeOld * 0.001 * 0;
% % Error 
%     for errorCount = 1:N-2
%         perturbVector = perturbVec(errorCount,1:2);
%         perturbVector = perturbVector/norm(perturbVector);
%         out(errorCount+1,1:2) = out(errorCount+1,1:2) + errorFactor * perturbVector;
%     end
%     
%     [out,~]=mvsplint(out,N);
% % 
% % %% Making drop axisymmetric
% % 
% out = errorReduction(out,N);
% [out,~] = mvsplint(out,N);
% Fitting spline over the digitised profile

outOld = mvsplint(pp,length(pp));
out = outOld;

% % Calculating original volume using pappus centroid theorem
volumeOld = pappus(out);
ReqOld = (volumeOld*3/(4*pi))^(1/3);
ReqNew = (3/(4*pi))^(1/3);
noramalisingFactor = ReqNew/ReqOld;

out_first(:,1:2) = out(:,1:2) * noramalisingFactor;
    
out = errorReduction(out_first);    % Makes drop profile symmetric about y axis

% % Smoothening drop profile 
% % Stepwise decrement of no. of coordinates in digitized drop profile to N
% NN = length(out);
% while NN > N
% 	out = mvsplint(out,NN);
% 	NN = NN - 1;
% end
% out = errorReduction(out);

out = mvsplint(out,N);
out = errorReduction(out);
out = mvsplint(out,N);

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

% Calculating original volume using pappus centroid theorem
volume_first = pappus(out);
Req = (volume_first*3/(4*pi))^(1/3);

%%
% Pre allocate matrices
meanCurvature = zeros(M1,N);
eqvCurvature = zeros(M1,N);
stdCurvature = zeros(M1,1);

if changeBo == 1
    % Considering droplet is under the effect of surface tension and gravitational force
    BoArray = zeros(M1,1);
    for count = 1:M1
        constant = Bo/(2*(Req^2));
        centrifugal = (1/4)*(Bw/Req^3);
        meanCurvature(count,1:N) = out(:,4);    % mean curvature of droplet
        eqvCurvature(count,1:N) = out(:,4) - out(:,2).*constant + (centrifugal)*(out(:,1).^2);    % equivalent curvature of droplet
        stdCurvature(count) = std(out(:,4) + (centrifugal)*(out(:,1).^2) - out(:,2).*constant);   % standard deviation of equivalent curvature
   
        BoArray(count) = Bo;
        Bo = Bo + BoIncrement;
    end
    
    % Bo value corresponding to minimum value of standard deviation in
    % equivalent curvature
    [minStd,minStdInd] = min(stdCurvature);
    Bomin = BoArray(minStdInd);
    
    % plot standard deviation vs guess values of Bo
    plot(BoArray(1:count),stdCurvature(1:count), 'k');
end

if changeBw == 1
% Considering droplet is under the effect of surface tension and centrifugal force
    BwArray = zeros(M2,1);
    for count = 1:M2
        constant = Bo/(2*(Req^2));
        centrifugal = (1/4)*(Bw/Req^3);
        meanCurvature(count,1:N) = out(:,4);    % mean curvature of droplet
        eqvCurvature(count,1:N) = out(:,4) - out(:,2).*constant + (centrifugal)*(out(:,1).^2);    % equivalent curvature of droplet
        stdCurvature(count) = std(out(:,4) + (centrifugal)*(out(:,1).^2) - out(:,2).*constant);   % standard deviation of equivalent curvature
   
        BwArray(count) = Bw;
        Bw = Bw + BwIncrement;        
    end
    
    % plot standard deviation vs guess values of Bw
    BwArray = transpose(BwArray);
    stdCurvature = transpose(stdCurvature);
    plot(BwArray(1:count),stdCurvature(1:count));
end

if changeBoBw == 1
% Considering droplet is under the effect of surface tension, gravitational and centrifugal force
    BwArray = zeros(M2,1);
    BoArray = zeros(M1,1);
    stdCurvatureMin = zeros(M2);
    for count1 = 1:M2
        for count = 1:M1
  
            constant = Bo/(2*(Req^2));
            centrifugal = (1/4)*(Bw/Req^3);
            meanCurvature(count,1:N) = out(:,4);    % mean curvature of droplet
            eqvCurvature(count,1:N) = out(:,4) - out(:,2).*constant + (centrifugal)*(out(:,1).^2);    % equivalent curvature of droplet
            stdCurvature(count) = std(out(:,4) + (centrifugal)*(out(:,1).^2) - out(:,2).*constant);   % standard deviation of equivalent curvature
   
            BoArray(count) = Bo;
            Bo = Bo + BoIncrement;

        end
        
        Bo = BoStart;
        BwArray(count1) = Bw;
        Bw = Bw + BwIncrement;
        stdCurvatureMin(count1) = min(stdCurvature);
    end
    
    BoArray = transpose(BoArray);
    stdCurvature = transpose(stdCurvature);
    BwArray = transpose(BwArray);
    stdCurvatureMin = transpose(stdCurvatureMin); 
end
