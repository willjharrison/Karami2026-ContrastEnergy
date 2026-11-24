function csfFilter = createCSF(x,fmax,vmax,beta,trunc)

csfFilter = 10.^(log10(vmax) - log10(2) * ((log10(x) - log10(fmax))/(log10(2*beta)/2)).^2);

csfFilter(log10(x) < log10(fmax) & csfFilter<((vmax)-(trunc))) = (vmax)-(trunc);
