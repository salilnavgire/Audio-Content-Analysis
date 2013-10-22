import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
from pylab import *

#function to generate signal in part I
def wave():
	f = np.logspace(2.69897004, 3.698970004, num=44100)
	x3 = cumsum(ones(44100))/44100
	print f
	print 44100/f
	print sum(cumsum(44100/f))
	y2 = (0.95)*(np.sin(np.cumsum((2*np.pi)*(f/(44100)))))
	#pl.xlim(0,1)
	pl.plot (x3,y2)

	plt.show()

#function to plot the spectrogram
def spectrogram(win, NFFT, hop, pad):
	
	Fs = 44.1*1000
	f = np.logspace(2.69897004, 3.698970004, num=44100)
	y2 = (0.95)*(np.sin(np.cumsum((2*np.pi)*(f/(44.1*1000)))))
	Pxx, freqs, bins, im = specgram(y2, NFFT=NFFT, Fs=Fs, window=win, noverlap=NFFT-hop, pad_to = pad)
	print freqs
	pl.xlim(0,1)
	pl.ylim(0,22500)
	pl.plot(Pxx,freqs)
	plt.show()

Fs = 44.1*1000
f = np.logspace(2.69897004, 3.698970004, num=44100)
y2 = (0.95)*(np.sin(np.cumsum((2*np.pi)*(f/(44.1*1000)))))


#PART I : generate a singnal 
wave()


#PART II : Plot spectrogram

# plot rectangular window
#spectrogram(np.ones(window length), NFFT, hop size, zero padding size)
spectrogram(np.ones(1024), 1024, 32, 1024)

#plot for hamming window
#spectrogram(window, NFFT, hop size, zero padding size)
#spectrogram(np.ones(256), 256, 128, 256)
#spectrogram(np.ones(256), 256, 32, 256)

spectrogram(np.hamming(256), 256, 256, 256)

#plot for blackmann window
#spectrogram(window,NFFT, hop size, zero padding size)
spectrogram(np.blackman(128), 128, 64, 128)

print 'end'


