function [gaussian, x] = calculateGaussian(mean, stdev, x)

%x=linspace(mean-stdev*3, mean+stdev*3, steps);
%x=-1:0.01:1;
gaussian=(1/(sqrt(2*pi)*stdev)).*exp(0-((mean-x).^2./(2*stdev^2)));

end

