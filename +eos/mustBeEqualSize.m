function mustBeEqualSize(a,b)
%MUSTBEEQUALSIZE Custom validation function
%   Test for equal size
if ~isequal(size(a),size(b))
    eid = 'Size:notEqual';
    msg = 'Size of first input must equal size of second input.';
    throwAsCaller(MException(eid,msg))
end
end

