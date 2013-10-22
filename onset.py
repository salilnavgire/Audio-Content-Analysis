import numpy as np
import matplotlib.pyplot as plt
from pylab import*
import wave
import struct
import sys
from scikits.audiolab import *
import marlib
import random
from datetime import datetime
import operator
import scipy
import numpy as np
import marlib.audiofile as AF
from marlib import eps
import math
import scipy.io


np.seterr(invalid='ignore')


def max_peak_pick(x,L):
    x = np.asarray(x)
    if x.ndim==1:
        x = x[:,np.newaxis]
    
    X = np.zeros(x.shape)
    pad = np.zeros(L)
    for r in range(x.shape[1]):
        y = np.concatenate([pad,x[:,r],pad],axis=0)
        res = np.zeros(y.shape)
        for n in range(L, len(y) - 2*L):
            z = y[n-L:n+L]
            if (z[L]>z[:L]).all() and (z[L]>=z[L:]).all():  
                res[n] = z[L]
        
        X[:,r] = res[L:-L]
    
    return X

def stft(x, fs, framesz, hop):
    framesamp = int(framesz*fs)
    hopsamp = int(hop*fs)
    w = scipy.hamming(framesamp)
    X = scipy.array([scipy.fft(w*x[i:i+framesamp]) 
                     for i in range(0, len(x)-framesamp, hopsamp)])
    return X

def spectrogram(x,framesize,hopsize,window=None,nfft=None):
    x = np.asarray(x)
    
    # Input Dim Check
    if x.ndim==2:
        x = x.mean(axis=1)
    
    # Window Sanity Check
    if window is None:
        window = np.ones(framesize)
    elif window.shape[0] != framesize:
        raise ValueError("Window length must equal %d (nfft)."%nfft)
    
    # NFFT check
    if not nfft is None:
        if nfft < framesize:
            raise ValueError("nfft must be >= framesize")
    
    X = [np.fft.rfft(window*x_m,nfft) for x_m in shingle(x,framesize,hopsize)]
    return 20.0*np.log10(np.abs(np.asarray(X)) + eps) 

def shingle(x,N,hopsize):
    M = int(np.ceil(x.shape[0] / float(hopsize)))
    X = np.zeros([M,N] + list(x.shape[1:]),dtype=x.dtype)
    for m in range(M):
        idx = int(np.round(hopsize*m))
        x_m = x[idx:idx+N]
        X[m,:x_m.shape[0]] = x_m
    
    return X



def princarg(ph):
	return np.pi + (ph + np.pi) % (-2*np.pi)


def generate_signal(filename):
	spf = wave.open(filename,'r')
	sound_info = spf.readframes(-1)
	sound_info = fromstring(sound_info, 'int32')
	p = spf.getnchannels()
	#print p
	f = spf.getframerate()
	#print f
	return sound_info, spf

def complex1(ph,Pxx):
	buff = shingle(ph, 128, 128)

	for i in range(0,len(buff)-1):
		for l in range(0,len(buff[0])-1):
			#print buff[i-1,j]
			#lol = 2*ph[i-1,j]-ph[i-2,j]
			ang = princarg(2*buff[i-1,l]-buff[i-2,l])*complex(0,1)
			Xk[i,l] = (Pxx[i-1,l])*np.exp(ang)
	return Xk


def Complex_Domain(Xk, Pxx):
	for i in range(0,len(Xk)-1):
		for j in range(0,len(Xk[0])-1):
			if (Pxx[i,j] - Xk[i,j]) < 0:
				CD[i] = 0
			else:
				CD[i] += np.absolute(Pxx[i,j] - Xk[i,j])
			#print Pxx[i,j] - Xk[i,j]
		CD[i] = CD[i]/128
	return CD


def peak(sal):
	for index, item in enumerate(sal):
	    if item > 0:
	        return index, item

def log_energy_derivative(sound_info):
	new = np.zeros(len(sound_info),dtype=float)
	new1 = np.zeros(len(sound_info),dtype=float)

	print len(sound_info)
	sqr = square(sound_info)
	pad = np.zeros(32)
	
	sound_pad = np.concatenate((pad,sqr,pad),axis=1)
	
	print "length of sound_pad is" + str(len(sound_pad))
	
	for i in range(0+len(pad),len(sound_info)):
		new[i-len(pad)] = (sum(np.absolute(sound_pad[i-(len(pad)):i+(len(pad))])))/len(pad)
		

	for i in range(1,len(new)-1):
		if np.isnan(new[i]-new[i-1]):
			new1[i] = new[i]
		if new[i]-new[i-1] == 0:
			new1[i] = 0
		else:
			new1[i] = (new[i]-new[i-1])

	npad = np.zeros(16)
	
	for i in range(0,len(new1)-len(npad)):
		new1[i] = sum(new[i:i+len(npad)])/len(npad)*np.average(new[i:i+len(npad)])

	
	"""
	new = np.concatenate((pad,new,pad),axis=1)
	for i in range(0,len(sound_info)-1):
		if new[i+len(pad)]-new[i+len(pad)-1]==0:
			new1[i] = 0
		else:
			new1[i] = np.log10(new[i+len(pad)]-new[i+len(pad)-1])

	"""
	
	return new1
	


def rectified_spectral_flux(sound_info,spf):

	f = spf.getframerate()
	Pxx, freqs, bins, im = specgram(sound_info, Fs = f, scale_by_freq=True,sides='default')
	pad = np.zeros(16)

	hfc = np.zeros(len(Pxx[0]),dtype=float)
	sf = np.zeros(len(Pxx[0]),dtype=float)

	Pxx.astype(float)
	row,col = shape(Pxx)
	print "shape is"
	print len(Pxx[0])

	for i in range(0,col-1):
		hfc[i] = sum(np.power(np.absolute(Pxx[:,i]),1)*i)/(col/2)
		
	#spectral flux
	for i in range(0,col-1):
		sf[i] = (sum(np.power(np.absolute((Pxx[:,i])-np.absolute(Pxx[:,i-1])+np.absolute(np.absolute(Pxx[:,i])-np.absolute(Pxx[:,i-1])))/2,1))/col)*2

	print max(sf)

	npad = np.zeros(4)

	for i in range(0,len(sf)-len(npad)-1):
		sf[i] = sum(sf[i:i+len(npad)])/len(npad)*np.average(sf[i:i+len(npad)])
	
	return sf


def peak_picking(func,win):
	sal = max_peak_pick(func, win)

	index = np.nonzero(sal)
	
	goli = index[0]


	amp = np.zeros(len(goli))
	for i in range(0,len(goli)):
		amp[i] = np.absolute(func[goli[i]])

	return index[0], amp

def vertcal_lines(func):

	mat = scipy.io.loadmat('tabla.mat')

	Time = spf.getnframes()/spf.getframerate()

	g = mat['T']

	#for tabla
	amp = np.zeros(len(g))
	
	gaa = np.zeros(len(g))
	
	#for piano
	"""
	amp = np.zeros(len(g[0]))
	
	gaa = np.zeros(len(g[0]))
	"""

	#for tabla
	for i in range(0,len(g)):
		#print g[i][0]
		gaa[i] = g[i][0]


	#for piano
	"""
	for i in range(len(g[0])):
		gaa[i] = g[0][i]
	"""

	for i in range(0,len(gaa)):
		gaa[i] = gaa[i]*(len(func)/Time)

	return gaa, amp

def stat(gaa, det):
	
	TP = 0
	FN = 0
	FP =0
	corr = 0


	for i in range(len(gaa)):

		for k in range(len(det)):
			if np.absolute(gaa[i]-det[k]) < 0.05:
				TP+=1
				break
			else:
				pass
		if k == (len(det)-1):
			FN+=1

	for i in range(len(det)):
		for k in range(len(gaa)):
			if np.absolute(det[i]-gaa[k]) < 0.05:
				corr+=1
				break
			else:
				pass
		if k == (len(gaa)-1):
			FP+=1


	print TP
	print FN
	print corr
	print FP


	P = TP/float(TP+FP)
	#print P
	R = TP/float(TP+FN)
	#print R
	F = 2*P*R/float(P+R+1)
	#print F

	return P, R, F
	


def convert(func,Time,index):
	sm = np.zeros(len(index))
	sm = sm.astype(np.float)

	for i in range(len(index)):
		
		sm[i] = float(index[i]/float(len(func)))*Time

	return sm	





"""Main Program"""

"""Part 1"""


print "Part 1"

"""Select Signal"""

spf = wave.open("tabla.wav",'r')
#print spf.getsampwidth()
sao = spf.readframes(-1)
sao = fromstring(sao, 'int16')
p = spf.getnchannels()
#print p

f = spf.getframerate()


sound_info = np.zeros(len(sao),dtype=float)

sao = sao.astype(np.float)

sound_info = sao/max(sao)

#print max(sound_info)

ph = np.zeros(len(sound_info),dtype=float)

Pxx = spectrogram(sound_info, 256, 128)

ph = (np.angle(sound_info))
print shape
buff = shingle(ph, 128, 128)
print shape(buff)


"""Log Energy derivative"""

new1 = np.zeros(len(sound_info),dtype=float)
new1 = log_energy_derivative(sound_info)


"""rectified spetral"""

sf = rectified_spectral_flux(sound_info,spf)


"""complex domain"""

Xk = np.empty((len(buff),len(buff[0])), dtype=np.complex)
ph = np.zeros(len(sound_info),dtype=float)
Pxx = spectrogram(sound_info, 256, 128)

ph = (np.angle(sound_info))

Xk = complex1(ph,Pxx)

CD = np.zeros(len(buff),dtype=float)

CD = Complex_Domain(Xk,Pxx)


#plots for part 1

"""
subplot(411)
plot(sound_info)
title('Signal')
index1, amp1 = peak_picking(sound_info,3000)
plot(index1, amp1, 'ko')
xlim(0,len(sound_info))


subplot(412)
plot(new1)
title('Log Energy Derivative')
index2, amp2 = peak_picking(new1,3000)
plot(index2, amp2, 'ko')
xlim(0,len(new1))


subplot(413)
plot(sf)
title('Rectified Spectral Flux')
index3, amp3 = peak_picking(sf,16)
plot(index3, amp3, 'ko')
xlim(0,len(sf))

subplot(414)
plot(CD)
title('Rectified Complex Domain')
index4, amp4 = peak_picking(CD,32)
plot(index4, amp4, 'ko')
xlim(0,len(CD))

"""

"""Part 2"""

print 'part2'

mat = scipy.io.loadmat('tabla.mat')
print "mat"
#print mat

Time = spf.getnframes()/spf.getframerate()

fr = spf.getnframes()

g = mat['T']

print g

#for tabla
amp = np.zeros(len(g))
	
gaa = np.zeros(len(g))


#for piano
"""
amp = np.zeros(len(g[0]))
	
gaa = np.zeros(len(g[0]))
"""


#for tabla
for i in range(0,len(g)):
	print g[i][0]
	gaa[i] = g[i][0]


#for piano
"""
for i in range(0,len(g[0])):
	print g[0][i]
	gaa[i] = g[0][i]
"""


print gaa
print len(gaa)

subplot(411)
plot(sound_info)
title('Signal')
index1, amp1 = peak_picking(sound_info,3000)
plot(index1, amp1, 'ko')
det = convert(sound_info,Time,index1)
P, R, F = stat(gaa,det)
print "Signal"
print "Precision "+str(P)
print "Recall "+str(R)
print "F-measure "+str(F)
g,amp = vertcal_lines(sound_info)
vlines(g,0,max(sound_info),'g')
xlim(0,len(sound_info))


subplot(412)
plot(new1)
title('Log Energy Derivative')
index2, amp2 = peak_picking(new1,3000)
plot(index2, amp2, 'ko')
det = convert(new1,Time,index2)
P, R, F = stat(gaa,det)
print "Log Energy Derivative"
print "Precision "+str(P)
print "Recall "+str(R)
print "F-measure "+str(F)
g,amp = vertcal_lines(new1)
vlines(g,0,max(new1),'g')
xlim(0,len(new1))


subplot(413)
plot(sf)
title('Rectified Spectral Flux')
index3, amp3 = peak_picking(sf,16)
plot(index3, amp3, 'ko')
det = convert(sf,Time,index3)
P, R, F = stat(gaa,det)
print "Rectified SPectral Flux"
print "Precision "+str(P)
print "Recall "+str(R)
print "F-measure "+str(F)
g,amp = vertcal_lines(sf)
vlines(g,0,max(sf),'g')
xlim(0,len(sf))

subplot(414)
plot(CD)
title('Rectified Complex Domain')
index4, amp4 = peak_picking(CD,32)
plot(index4, amp4, 'ko')
det = convert(CD,Time,index4)
P, R, F = stat(gaa,det)
print "Rectified Complex Domain"
print "Precision "+str(P)
print "Recall "+str(R)
print "F-measure "+str(F)
g,amp = vertcal_lines(CD)
vlines(g,0,max(CD),'g')
xlim(0,len(CD))

show()


print "end"
