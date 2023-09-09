function mustBeCorrectionFactor(arg)
%MUSTBECORRECTIONFACTOR Checks if arg is a correction factor
if ~isa(arg,'eos.UnitCorrectionFactor') && ~isa(arg,'eos.SoaveCorrectionFactor')
    eid = 'Type:notCorrectionFactor';
    msg = 'Type of arg must be a type of correction factors.';
    throwAsCaller(MException(eid,msg))
end
end

