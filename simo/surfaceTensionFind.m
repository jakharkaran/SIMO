function [out, surfaceTension, volumeScaled, ReqScaled  ] = surfaceTensionFind( pp, AlgorithmBo, density, scalingQuantity, scale, scalingQuantityPer )
%% Surface Tension Calculator
% Determine surface tension using Bo and length scale (pixelpermm)

N = 50;
Bo = AlgorithmBo;
g = 9.81;

outOld = mvsplint(pp,length(pp));
out = outOld;

out = errorReduction(out);
NN = length(out);
while NN > N
	out = mvsplint(out,NN);
	NN = NN - 1;   
end

out = errorReduction(out);
outNew = mvsplint(out,N);
out = outNew;
volumeUnScaled = pappus(outNew);
ReqUnScaled = (volumeUnScaled*3/(4*pi))^(1/3); 

if strcmp(scale, 'cubic metre') || strcmp(scale, 'microlitre') || strcmp(scale, 'nanolitre')
    
    if strcmp(scale, 'cubic metre')
        volumeScaled = scalingQuantity;
    elseif strcmp(scale, 'microlitre')
        volumeScaled = scalingQuantity * 10^-9;
    elseif strcmp(scale, 'nanolitre')
        volumeScaled = scalingQuantity * 10^-12;
    end
    
    ReqScaled = (volumeScaled*3/(4*pi))^(1/3);   
    normalizingFactor = ReqScaled/ReqUnScaled;   
    out(:,1:2) = out(:,1:2) * normalizingFactor;
    out = errorReduction(out);
    out = mvsplint(out,N);
    
elseif strcmp(scale, 'pixel per milimetre') || strcmp(scale, 'pixel per metre')

    if strcmp(scale, 'pixel per milimetre') || strcmp(scale, 'milimetre')
       out = out/(scalingQuantity * 10^3);       
    elseif strcmp(scale, 'pixel per metre') || strcmp(scale, 'metre')
       out = out/scalingQuantity;
    end
    
    out = errorReduction(out);
    out = mvsplint(out,N);
    volumeScaled = pappus(out);
    ReqScaled = (volumeScaled*3/(4*pi))^(1/3);
    
elseif  strcmp(scale, 'milimetre') || strcmp(scale, 'metre')
    
    if strcmp(scale, 'milimetre')
        out = out/(scalingQuantity * 10^3 / scalingQuantityPer); 
    elseif strcmp(scale, 'metre')
        out = out/(scalingQuantity / scalingQuantityPer); 
    end
    
    out = errorReduction(out);
    out = mvsplint(out,N);
    volumeScaled = pappus(out);
    ReqScaled = (volumeScaled*3/(4*pi))^(1/3);
end

surfaceTension = density*g*ReqScaled^2/Bo;

end


