function isoct = isoctave()
%Michael Hirsch
% tested with Octave 3.6-4.2 and Matlab

persistent oct;

if isempty(oct)
    oct = exist('OCTAVE_VERSION', 'builtin') == 5;
end

isoct=oct; % has to be a separate line/variable for matlab

end
