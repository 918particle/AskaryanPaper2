function res = fit_function(p0,x)
	xfit = (x-p0(1))./p0(2);
	res = p0(3).*xfit.^2.*exp(-0.5*xfit.^2);
endfunction