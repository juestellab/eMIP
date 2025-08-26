function [alphac, Pj] = FBE_alpha(LF, HF, opt)
%FBE_ALPHA: Frequency Band Equalizer (see: DOI: 10.1038/s41551-017-0068)
    arguments
        LF                                                                                                                      = [];
        HF                                                                                                                      = [];
        opt.FBEoptions     char        {mustBeMember(opt.FBEoptions, {'Pjx', 'Pjy', 'Pjz', 'Pji', 'volume','energy'})}        = 'Pjz';                        % Pjx_y_z: FBE based on Pj along x,y,z axis; ; Pji: FBE for each Pj individually; volume: FBE based on 3D volume; energy; FBE based on volume squared  
        opt.Pj              struct                                                                                              = struct('LF',[],'HF',[]);      % predefined Pjs (if available) to avoid unnecessary projections   
        opt.PjType          char        {mustBeMember(opt.PjType, {'MIP', 'Energy', 'Sum'})}                                    = 'MIP';
    end
    
    switch lower(opt.FBEoptions)
    case {'pjz', 'pjx', 'pjy'}
        pjDir      = lower(opt.FBEoptions(3));
        axisDir     = pjDir2axisDir(pjDir);

        if ~isfield(opt.Pj.LF, pjDir)
            opt.Pj.LF.(pjDir)      = Vol2ImPj(LF, "AxisDir", axisDir, "ProjectionType", opt.PjType);  
        end
        LF          = opt.Pj.LF.(pjDir);

        if ~isfield(opt.Pj.HF, pjDir)
            opt.Pj.HF.(pjDir)      = Vol2ImPj(HF, "AxisDir", axisDir, "ProjectionType", opt.PjType);  
        end
        HF          = opt.Pj.HF.(pjDir);
    case 'energy'
        HF          = HF.^2;
        LF          = LF.^2;
    case 'pji'
        alphac = 1;
        return
    end 

    alphac = lsqlin(double(HF(:)),double(LF(:)));
    
    switch lower(opt.FBEoptions)
    case 'energy'
        alphac      = sqrt(alphac);
    end
    Pj          = opt.Pj;
end