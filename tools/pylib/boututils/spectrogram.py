"""
Creates spectrograms using the Gabor transform to maintain time and frequency resolution

(spectrogram, frequency) = spectrogram(data, dt, sigma)

four arguments are required:

	data  - the time series you want spectrogrammed
	dt    - time resolution
	sigma - used in the Gabor transform, will balance time and frequency resolution
		suggested value is 1.0, but may need to be adjusted manually
		until result is as desired
		
		IF bands are too tall raise sigma
		IF bands are too wide, lower sigma
		
optional keyword:

	clip = 1.0  - optional keyword, makes the spectrogram run faster, but decreases frequency resolution
		clip is by what factor the time spectrum should be clipped by --> N_new = N / clip
	
	plot = True - do you want it to plot the data for you? (True/False)
	
	
NOTE: Very early and very late times will have some issues due to the method - truncate them after
	taking the spectrogram if they are below your required standards
	
NOTE: If you are seeing issues at the top or bottom of the frequency range, you need a longer time series

written by: Jarrod Leddy
updated:    4/10/2013

"""

try:
	from numpy import arange,zeros,exp,power,transpose,sin,cos,linspace,min,max
	from scipy import fftpack,pi
	from scipy import pi
	import matplotlib.pyplot as plt

except ImportError:
	print "ERROR: NumPy or SciPy module not available"
	raise
	
def spectrogram(data, dx, sigma, clip=1.0, plot=True):
	n = data.size
	xx = arange(n)*dx
	sigma = sigma * 0.01 * dx
	
	if clip: # clip to 6 std dev on either side for speed
		
		n_clipped = int(n/clip)
		
		spectra = zeros((n,n_clipped/2))
	
		omega = fftpack.fftfreq(n_clipped, dx)
		omega = omega[0:n_clipped/2]
	
		for i in range(n):
			beg = i-n_clipped/2
			end = i+n_clipped/2-1
			if beg < 0:
				end = end-beg
				beg = 0
			elif end >= n:
				end = n-1
				beg = end - n_clipped + 1 
			gaussian = 1.0 / (sigma * 2.0 * pi) * exp(-0.5 * power(( xx[beg:end] - xx[i] ),2.0) / (2.0 * sigma) )
			fftt = abs(fftpack.fft(data[beg:end] * gaussian))
			fftt = fftt[:n_clipped/2]
			spectra[i,:] = fftt
		
	else:
	
		spectra = zeros((n,n/2))
		
		omega = fftpack.fftfreq(n, dx)
		omega = omega[0:n/2]
	
		for i in range(n):
			gaussian = 1.0 / (sigma * 2.0 * pi) * exp(-0.5 * power(( xx - xx[i] ),2.0) / (2.0 * sigma) )
			fftt = abs(fftpack.fft(data * gaussian))
			fftt = fftt[:n/2]
			spectra[i,:] = fftt
	
	if plot:
		try:
			import matplotlib.pyplot as plt
		except ImportError:
			print "ERROR: MatPlotLib module not available"
			raise
		
		
		nt = len(data)
		print arange(nt).shape, omega.shape, spectra.shape
		lvls = linspace(min(spectra),max(spectra),100)
		plt.contourf(arange(nt)*dx,omega,spectra.T,levels=lvls)
		plt.ylabel("frequency")
		plt.xlabel(r"$t$")
		plt.title(r"Spectrogram")
		plt.show()
	
	return (transpose(spectra), omega)
	
def test_spectrogram(n, d, s):
	
	"""
	Function used to test the performance of spectrogram with various values of sigma
	"""
	xx = arange(n)/d
	test_data = sin(2.0*pi*512.0*xx * ( 1.0 + 0.005*cos(xx*50.0))) + 0.5*exp(xx)*cos(2.0*pi*100.0*power(xx,2))
	test_sigma = s
	dx = 1.0/d
	
	s1 = test_sigma*0.1
	s2 = test_sigma
	s3 = test_sigma*10.0
	
	(spec2,omega2) = spectrogram(test_data, dx, s2, clip=1.0, plot=False)
	(spec3,omega3) = spectrogram(test_data, dx, s3, clip=1.0, plot=False)
	(spec1,omega1) = spectrogram(test_data, dx, s1, clip=1.0, plot=False)
	
	levels = linspace(min(spec1),max(spec1),100)
	plt.subplot(311)
	plt.contourf(xx,omega1,spec1,levels=levels)
	plt.ylabel("frequency")
	plt.xlabel(r"$t$")
	plt.title(r"Spectrogram of $sin(t + cos(t) )$ with $\sigma=$%3.1f"%s1)
	
	levels = linspace(min(spec2),max(spec2),100)
	plt.subplot(312)
	plt.contourf(xx,omega2,spec2,levels=levels)
	plt.ylabel("frequency")
	plt.xlabel(r"$t$")
	plt.title(r"Spectrogram of $sin(t + cos(t) )$ with $\sigma=$%3.1f"%s2)
	
	levels = linspace(min(spec3),max(spec3),100)
	plt.subplot(313)
	plt.contourf(xx,omega3,spec3,levels=levels)
	plt.ylabel("frequency")
	plt.xlabel(r"$t$")
	plt.title(r"Spectrogram of $sin(t + cos(t) )$ with $\sigma=$%3.1f"%s3)
	plt.show()
	
if __name__ == "__main__":
	test_spectrogram(2048, 2048.0, 1.0) # array size, divisions per unit, sigma of gaussian
