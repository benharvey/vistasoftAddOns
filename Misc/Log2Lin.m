fwhms=model{1}.sigma.major.*(2*sqrt(2*log(2)));
fwhms=exp(model{1}.x0+fwhms./2)-exp(model{1}.x0-fwhms./2);
model{1}.sigma.major=fwhms;
model{1}.sigma.minor=fwhms;


if isfield(model{1}, 'sigma2')
    fwhms2=model{1}.sigma2.major.*(2*sqrt(2*log(2)));
    fwhms2=exp(model{1}.x0+fwhms2./2)-exp(model{1}.x0-fwhms2./2);
    model{1}.sigma2.major=fwhms2;
    model{1}.sigma2.minor=fwhms2;
end

model{1}.x0=exp(model{1}.x0);