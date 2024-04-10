clear;

Bo = 0.2;
BoStart = -10;
BoEnd = 10  ;

count = 1;
counter = 1;

N = 20;

BoAssume = [BoStart ((BoStart + ((BoStart + BoEnd)/2))/2) (BoStart + BoEnd)/2 ((BoEnd + ((BoStart + BoEnd)/2))/2) BoEnd];

while(BoAssume(5) - BoAssume(1) > 0.001)
% for count = 1:N
        
    BoAssumeOld = BoAssume;
    
    X = rand;
    if (X < (1/3))
%     if  (Bo >= BoAssume(1)) && (Bo < ((abs(BoAssume(3) - BoAssume(2)) * (1/3)) + BoAssume(2)))
        BoAssume(1) = BoAssumeOld(1);
        BoAssume(3) = BoAssumeOld(2);
        BoAssume(5) = BoAssumeOld(3);
        BoAssume(2) = (BoAssume(1)+BoAssume(3))/2;
        BoAssume(4) = (BoAssume(3)+BoAssume(5))/2;

    elseif (X>=(1/3)) && (X<=(1/3))
%     elseif (Bo >= ((abs(BoAssume(3) - BoAssume(2)) * (1/3)) + BoAssume(2))) && (Bo <= ((abs(BoAssume(4) - BoAssume(3)) * (2/3)) + BoAssume(3)))
        BoAssume(1) = BoAssumeOld(2);
        BoAssume(3) = BoAssumeOld(3);
        BoAssume(5) = BoAssumeOld(4);
        BoAssume(2) = (BoAssume(1)+BoAssume(3))/2;
        BoAssume(4) = (BoAssume(3)+BoAssume(5))/2;

    elseif (X>(1/3))
%     elseif (Bo > ((abs(BoAssume(4) - BoAssume(3)) * (2/3)) + BoAssume(3))) && (Bo <= BoAssume(5))
        BoAssume(1) = BoAssumeOld(3);
        BoAssume(3) = BoAssumeOld(4);
        BoAssume(5) = BoAssumeOld(5);
        BoAssume(2) = (BoAssume(1)+BoAssume(3))/2;
        BoAssume(4) = (BoAssume(3)+BoAssume(5))/2;

    end
Bo5(count, :) = BoAssumeOld;
count = count + 1;
end

BoAssume = [BoStart (BoStart+(abs(BoStart-BoEnd)*(1/3))) (BoEnd-(abs(BoStart-BoEnd)*(1/3))) BoEnd];

%% Golden Swction Method
% for counter = 1:N

while((BoAssume(4) - BoAssume(1)) > 0.001)
    BoAssumeOld = BoAssume;
    
    X = rand;
    
    if X < 0.5
%     if (Bo >= BoAssume(1)) && (Bo < ((BoAssume(2) + BoAssume(3))/2))
        BoAssume(1) = BoAssumeOld(1);
        BoAssume(4) = BoAssumeOld(3);
        BoAssume(2) = ((abs(BoAssume(1)-BoAssume(4))/3) + BoAssume(1));
        BoAssume(3) = (BoAssume(4) - (abs(BoAssume(1)-BoAssume(4))/3));
        
    elseif X>=0.5
%     elseif (Bo <= BoAssume(4)) && (Bo >= ((BoAssume(2) + BoAssume(3))/2))
        BoAssume(1) = BoAssumeOld(2);
        BoAssume(4) = BoAssumeOld(4);
        BoAssume(2) = ((abs(BoAssume(1)-BoAssume(4))/3) + BoAssume(1));
        BoAssume(3) = (BoAssume(4) - (abs(BoAssume(1)-BoAssume(4))/3));       
%     end
    end
    
Bo4(counter,:) = BoAssume;
counter = counter + 1;
end

% h = figure;
% plot(1:N, Bo5(:,5) - Bo5(:,1), '.r');
% hold on;
% plot(1:N, Bo4(:,4) - Bo4(:,1), '.b');
% hold off;

axis([-1 N+1 -1 BoEnd-BoStart+1]);
