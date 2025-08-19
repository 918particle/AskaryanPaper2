function res = fit_function_3(p0,x)
	xfit = (x-p0(1))./p0(2);
	res = p0(3).*exp(-0.5*xfit.^2);
endfunction