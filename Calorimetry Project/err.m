%purpose: Calculates largest value of sigma y and returns to main function
%inputs: Q matrix from each fitline functino call, time of interest,
%calculated sigma y from fitline function
%data: NA
%outputs: final sigma_y
%assumptions: N/A
%Author names: Matt Bechtel, Pax Armon
%ID numbers: 109802403, 109799973
%Date created - 10/12/21
%Date modified - 10/14/21

function sigma_y_new = err(Q1, Q2, t, sigma_y_1)

%clean up matrix math
a = [t 1];
b = [t; 1];

%calculate sigma y associated with each Q matrix
sigma_y_21 = sqrt(a * Q1 * b);
sigma_y_22 = sqrt(a * Q2 * b);

%combine sigma y from each Q matrix in quadrature
sigma_y = sqrt((sigma_y_21)^2 + (sigma_y_21)^2);

%determine largest sigma y and assigns to variable that is returned to main
%function
if sigma_y > sigma_y_1
    sigma_y_new = sigma_y(1);
elseif sigma_y_1 > sigma_y
    sigma_y_new = sigma_y_1;
else
    sigma_y_new = sigma_y;
end
end