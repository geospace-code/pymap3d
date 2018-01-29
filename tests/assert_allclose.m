function ok = assert_allclose(actual, desired, atol, rtol)
% ok = assert_allclose(actual, desired, atol, rtol)
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

    if nargin < 4, rtol=1e-5; end
    if nargin < 3 || isempty(atol), atol=1e-8; end
    
    assert(isscalar(rtol))
    assert(isscalar(atol))

    [maxerrmag,i] = max(abs(actual-desired));
    if maxerrmag > atol + rtol * abs(desired)
      disp(['Actual: ',num2str(actual)])
      disp([' desired: ',num2str(desired)])
      error(['AssertionError: maximum error magnitude ',num2str(maxerrmag),' on ',int2str(i),'th value: ',num2str(actual(i)),' atol: ',num2str(atol),' rtol: ',num2str(rtol)])
    end

end

% Copyright (c) 2014-2018 Michael Hirsch, Ph.D.
%
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
