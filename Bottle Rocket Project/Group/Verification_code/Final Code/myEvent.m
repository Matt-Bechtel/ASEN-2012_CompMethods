% Purpose: This function sets the limit on the ODE45 integration. Learned
%          this in Office Hours. Not quite sure how it works.
% 
% Inputs: NA
%
% Outputs: NA
%
% Assumptions: NA
%
% Author: Matt Bechtel
%
% ID Number: 109802403
%
% Date Created: 11/9/21
%
% Date Modified: 11/19/21

function [val, ist, dur] = myEvent(~, X)
    val = (X(2) <= 0);
    ist = 1;
    dur = 0;
end