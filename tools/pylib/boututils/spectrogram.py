"""Creates spectrograms using the Gabor transform to maintain time and
frequency resolution

written by: Jarrod Leddy
updated:    23/06/2016

"""
from __future__ import print_function
from __future__ import division
from builtins import range

from numpy import arange, zeros, exp, power, transpose, sin, cos, linspace, min, max
from scipy import fftpack, pi


def spectrogram(data, dx, sigma, clip=1.0, optimise_clipping=True, nskip=1.0):
    """Creates spectrograms using the Gabor transform to maintain time
    and frequency resolution

    .. note:: Very early and very late times will have some issues due
          to the method - truncate them after taking the spectrogram
          if they are below your required standards

    .. note:: If you are seeing issues at the top or bottom of the
          frequency range, you need a longer time series

    written by: Jarrod Leddy
    updated:    23/06/2016

    Parameters
    ----------
    data : array_like
        The time series you want spectrogrammed
    dt : float
        Time resolution
    sigma : float
        Used in the Gabor transform, will balance time and frequency
        resolution suggested value is 1.0, but may need to be adjusted
        manually until result is as desired:

            - If bands are too tall raise sigma
            - If bands are too wide, lower sigma
    clip : float, optional
        Makes the spectrogram run faster, but decreases frequency
        resolution. clip is by what factor the time spectrum should be
        clipped by --> N_new = N / clip
    optimise_clip : bool
        If true (default) will change the data length to be 2^N
        (rounded down from your inputed clip value) to make FFT's fast
    nskip : float
        Scales final time axis, skipping points over which to centre
        the gaussian window for the FFTs

    Returns
    -------
    tuple : (array_like, array_like, array_like)
        A tuple containing the spectrogram, frequency and time

    """
    n = data.size
    nnew = int(n/nskip)
    xx = arange(n)*dx
    xxnew = arange(nnew)*dx*nskip
    sigma = sigma * dx

    n_clipped = int(n/clip)

    # check to see if n_clipped is near a 2^n factor for speed
    if(optimise_clipping):
        nn = n_clipped
        two_count = 1
        while(1):
            nn = nn/2.0
            if(nn <= 2.0):
                n_clipped = 2**two_count
                print('clipping window length from ',n,' to ',n_clipped,' points')
                break
            else:
                two_count += 1
    else:
        print('using full window length: ',n_clipped,' points')

    halfclip = int(n_clipped/2)
    spectra = zeros((nnew,halfclip))

    omega = fftpack.fftfreq(n_clipped, dx)
    omega = omega[0:halfclip]

    for i in range(nnew):
        beg = i*nskip-halfclip
        end = i*nskip+halfclip-1

        if beg < 0:
            end = end-beg
            beg = 0
        elif end >= n:
            end = n-1
            beg = end - n_clipped + 1

        gaussian = 1.0 / (sigma * 2.0 * pi) * exp(-0.5 * power(( xx[beg:end] - xx[i*nskip] ),2.0) / (2.0 * sigma) )
        fftt = abs(fftpack.fft(data[beg:end] * gaussian))
        fftt = fftt[:halfclip]
        spectra[i,:] = fftt

    return (transpose(spectra), omega, xxnew)


def test_spectrogram(n, d, s):
    """Function used to test the performance of spectrogram with various
    values of sigma

    Parameters
    ----------
    n : int
        Number of points
    d : float
        Grid spacing
    s : float
        Initial sigma

    """

    import matplotlib.pyplot as plt

    nskip = 10
    xx = arange(n)/d
    test_data = sin(2.0*pi*512.0*xx * ( 1.0 + 0.005*cos(xx*50.0))) + 0.5*exp(xx)*cos(2.0*pi*100.0*power(xx,2))
    test_sigma = s
    dx = 1.0/d

    s1 = test_sigma*0.1
    s2 = test_sigma
    s3 = test_sigma*10.0

    (spec2,omega2,xx) = spectrogram(test_data, dx, s2, clip=5.0, nskip=nskip)
    (spec3,omega3,xx) = spectrogram(test_data, dx, s3, clip=5.0, nskip=nskip)
    (spec1,omega1,xx) = spectrogram(test_data, dx, s1, clip=5.0, nskip=nskip)

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
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    test_spectrogram(2048, 2048.0, 0.01) # array size, divisions per unit, sigma of gaussian
