%purpose: This function takes a weighted average of two independent
%variables. It was used to take a weighted average of slopes and intercepts
%to combine LLS lines into one line.
%inputs: two variable to be averaged together, error associated with each
%variable
%data: NA
%outputs: weighted average, error associated with weighted average 
%assumptions: NA
%Author names: Matt Bechtel, Pax Armon
%ID numbers: 109802403, 109799973
%Date created - 10/5/21
%Date modified - 10/14/21

function [x_wav,sigma_wav] = weighted_average(x1 ,x2, sigma_x1, sigma_x2) % Weighted average function definition: take in a matrix
 % Calculate the combined weighted average of the data columns
 weight_x1 = 1 / sigma_x1^2;
 weight_x2 = 1 / sigma_x2^2; 
 
 % Calculate the uncertainty in the weighted average
 sigma_wav = 1 / sqrt((weight_x1 + weight_x2)); %good

 %combine each variable with its weight
 wav1 = x1 * weight_x1;
 wav2 = x2 * weight_x2;
 
 x_wav = (sum(wav1) + sum(wav2)) / (weight_x1 + weight_x2);
end