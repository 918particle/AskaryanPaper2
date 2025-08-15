function res = fit_function_2(p0,x)
	xfit = (x-p0(1))./p0(2);
	res = p0(2).*exp(-xfit);
endfunction