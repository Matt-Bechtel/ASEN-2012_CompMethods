%purpose: This function takes in data to create a system of matrices to
%determine the best estimate for slope and intercept for linear least
%squares regression. This also calculates one value for sigma y
%inputs: slope coefficients, solutions, time vector
%data: NA
%outputs: fit line, best estimate for m and b, error associated with m and
%b, sigma y, Q matrix
%assumptions: NA
%Author names: Matt Bechtel, Pax Armon
%ID numbers: 109802403, 109799973
%Date created - 10/4/21
%Date modified - 10/14/21

function [y_new, m_best, b_best, sigma_m_best, sigma_b_best, sigma_y_best, Q] = fitline(m_data, solutions, time)
    %% Linear Least Squares
    d = [solutions];

    %create column vectors for A matrix
    Acol1 = [m_data];
    Acol2 = ones(length(m_data), 1);

    %concatenate into one matrix
    A = [Acol1, Acol2];

    %solve for m and b
    x_hat = ([A' * A]^-1) * ( A' * d);
    
    %extract m and b from x hat matrix
    m_best = x_hat(1);
    b_best = x_hat(2);
   
    %new solution line
    y_new = m_best.*time + b_best;
    
    %% Error
    %clean up calculation
    a = 1 / (length(m_data) - 2);
    b = sum(((d - b_best - (m_best .* time)).^2));
    
    %calculate sigma y (book equation)
    sigma_y_best = sqrt(a * b);
    
    %create weight matrix of n x n dimensions that length of A matrix
    weight = eye(length(A)) * (1/sigma_y_best^2);
    
    %calculate covariance matrix
    Q = ([A' * weight * A]^-1);
    
    %find error associated with m and b
    sigma_m_best = sqrt(Q(1,1));
    sigma_b_best = sqrt(Q(2,2));
end