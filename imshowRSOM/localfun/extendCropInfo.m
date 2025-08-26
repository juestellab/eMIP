function [opt] = extendCropInfo(opt)
CropFields ={'Crop', 'CropInmm', 'CropInum', 'CropInPixel'};

for k=1:4
    if isempty(opt.(CropFields{k}))
        continue
    end
    j = 3-size(opt.(CropFields{k}),2);
    opt.(CropFields{k})(1:2,end+1:3) = repmat([0;inf],1,j);
end
end

