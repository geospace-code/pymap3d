function ok = assert_allclose(actual, desired, atol, rtol)
% atol: absolute tolerance
% rtol: relative tolerance
% based on numpy.testing.assert_allclose
% https://github.com/numpy/numpy/blob/v1.13.0/numpy/core/numeric.py#L2522
% for Matlab and GNU Octave
%
% if "actual" is within atol OR rtol of "desired", no error is emitted.

if nargin < 4, rtol=1e-5; end
if nargin < 3 || isempty(atol), atol=1e-8; end
  
maxerrmag = max(abs(actual-desired));
if maxerrmag > atol + rtol * abs(desired)
  error(['AssertionError: maximum error magnitude ',num2str(maxerrmag),' atol: ',num2str(atol),' rtol: ',num2str(rtol)])
end

end