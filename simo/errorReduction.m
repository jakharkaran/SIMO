function [ outNew ] = errorReduction( out )
%errorReduction Makes curve symmetric about y axis
    
N = length(out(:,1));
    
    if (rem(N,2)==0)
        % For even values of N
        
        for errorAverageCount = 1:N/2
            % Makes x axis coordinate symmetric about y axis
            meanxError = (out(errorAverageCount,1) + out(N+1-errorAverageCount,1))/2;
            out(errorAverageCount,1) = out(errorAverageCount,1) - meanxError;
            out(N+1-errorAverageCount,1) = out(N+1-errorAverageCount,1) - meanxError;
            meanXError(errorAverageCount) = meanxError;

            % Makes y axis coordinates symmetric about y axis
            meanyError = (out(errorAverageCount,2) + out(N+1-errorAverageCount,2))/2;
            out(errorAverageCount,2) = meanyError;
            out(N+1-errorAverageCount,2) = meanyError;
            meanYError(errorAverageCount) = meanyError;
        end

        sortedy = sort(out(1:2,2));
        out(1:2,2) = sortedy(1:2);
        out(N,2) = out(1,2);
        out(N-1,2) = out(2,2);
    %     out(N-2,3) = out(3,2);
    %     out(N-3,2) = out(4,2);
% 
%         for errorCount = 1:N/2
%             % Makes y coordinate of first and last contact point = 0.
%             out(errorCount,2) = out(errorCount,2) - out(1,2);
%             out(N+1-errorCount,2) = out(N+1-errorCount,2) - out(N,2);
%         end
%         
    else
        % For odd values of N
        
        for errorAverageCount = 1:(N-1)/2
            % Makes x axis coordinate symmetric about y axis
            meanxError = (out(errorAverageCount,1) + out(N+1-errorAverageCount,1))/2;
            out(errorAverageCount,1) = out(errorAverageCount,1) - meanxError;
            out(N+1-errorAverageCount,1) = out(N+1-errorAverageCount,1) - meanxError;
            meanXError(errorAverageCount) = meanxError;

            % Makes y axis coordinates symmetric about y axis
            meanyError = (out(errorAverageCount,2) + out(N+1-errorAverageCount,2))/2;
            out(errorAverageCount,2) = meanyError;
            out(N+1-errorAverageCount,2) = meanyError;
            meanYError(errorAverageCount) = meanyError;
        end
        out((N+1)/2,1) = 0;
            
        sortedy = sort(out(1:2,2));
        out(1:2,2) = sortedy(1:2);
        out(N,2) = out(1,2);
        out(N-1,2) = out(2,2);
    %     out(N-2,3) = out(3,2);
    %     out(N-3,2) = out(4,2);
    
    end
    
	% Makes y coordinate of first and last contact point = 0.
    out(:,2) = out(:,2) - out(1,2);
    outNew = out;
end

