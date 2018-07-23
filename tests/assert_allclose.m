function ok = assert_allclose(actual, desired, rtol, atol, equal_nan, err_msg,notclose,verbose)
% ok = assert_allclose(actual, desired, rtol, atol)
%
% Inputs
% ------
% atol: absolute tolerance (scalar)
% rtol: relative tolerance (scalar)
%
% Output
% ------
% ok: logical TRUE if "actual" is close enough to "desired"
%
% based on numpy.testing.assert_allclose
% https://github.com/numpy/numpy/blob/v1.13.0/numpy/core/numeric.py#L2522
% for Matlab and GNU Octave
%
% if "actual" is within atol OR rtol of "desired", no error is emitted.
narginchk(2,7)

if nargin < 8, verbose=false; end
if nargin < 7, notclose=false; end
if nargin < 6, err_msg=''; end
if nargin < 5 || isempty(equal_nan), equal_nan=false; end
if nargin < 4 || isempty(atol), atol=0; end
if nargin < 3 || isempty(rtol), rtol=1e-8; end

validateattributes(actual, {'numeric'}, {'real'})
validateattributes(desired, {'numeric'}, {'real'})
validateattributes(rtol, {'numeric'}, {'real','scalar','finite'})
validateattributes(atol, {'numeric'}, {'real','scalar','finite'})
validateattributes(equal_nan, {'numeric','logical'}, {'scalar'})
validateattributes(notclose, {'numeric','logical'}, {'scalar'})
validateattributes(verbose, {'numeric','logical'}, {'scalar'})
%% compute
  actual = actual(:);
  desired = desired(:);
  
  if equal_nan
    match = false(size(actual));
    match(isnan(actual)) = true;
  else
    match = true;
  end

  
  measdiff = abs(actual-desired);
  tol = atol + rtol * abs(desired);
  result = measdiff <= tol;
%% assert_allclose vs assert_not_allclose
  if notclose % more than N % of values should be changed more than tolerance (arbitrary)
    testok = match | sum(~result) > 0.0001*numel(desired);
  else
    testok = all(match | result);
  end
  
  if ~testok
    Nfail = sum(~result);
    j = find(~result);
    [bigbad,i] = max(measdiff(j));
    i = j(i);
    if verbose
      disp(['error mag.: ',num2str(measdiff(j)')])
      disp(['tolerance:  ',num2str(tol(j)')])
      disp(['Actual:     ',num2str(actual(i))])
      disp(['desired:    ',num2str(desired(i))])
    end

    error(['AssertionError: ',err_msg,' ',num2str(Nfail/numel(desired)*100,'%.2f'),'% failed accuracy. maximum error magnitude ',num2str(bigbad),' Actual: ',num2str(actual(i)),' Desired: ',num2str(desired(i)),' atol: ',num2str(atol),' rtol: ',num2str(rtol)])
  end

end

% Copyright (c) 2014-2018 Michael Hirsch, Ph.D.
%
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
