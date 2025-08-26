function [I] = visualizeRSOMFreqSep(img_LF,img_HF, dim, con_high, con_low)
%VISUALIZERSOMFREQSEP enables frequency seperated 2D visualization of RSOM
%   img_LF      low frequency component of RSOM image(3D)
%   img_HF      high frequency component of RSOM image(3D)
%   dim         dimension of projection for MIP
%   con_high    lower value for imadjust
%   con_low     higher value for imadjust

% check for valid input and set default values if required
switch nargin
    case 2
        dim = 2;
        con_high = 0.06;
        con_low = 0.35;
    case 3
        con_high = 0.06;
        con_low = 0.35;
    case 4
        con_low = 0.35;
    case 5
    otherwise
        error('Invalid number of input values')
end

I_low=squeeze(max(img_LF,[],dim));
I_high=squeeze(max(img_HF,[],dim));
[alphaval] = alphacalc(I_low,I_high);
I = imfuse(I_low,alphaval*I_high,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]); %[1 2 0]
I = imadjust(I,[con_high con_high 0; con_low con_low 1],[]); 
end

function [ alphac] = alphacalc(projlow,projhigh)

 alphaval=[0:0.01:100];
 differ=zeros(size(alphaval));
 for(i=1:length(alphaval))
 differ(i)=sum(sum((projlow - alphaval(i)*projhigh).^2));
 end
  [dum,ind]=min(differ);
  alphac=alphaval(ind);
%figure; plot(alphaval,differ);
end