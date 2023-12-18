from matplotlib import pyplot as plt
import pickle

#with open('objs.pkl','rb') as f:  # Python 3: open(..., 'rb')
f = open('objs.pkl','rb')
time0, doppler, t_end_of_pass, signal_freq = pickle.load(f)
f.close()
plt.figure(1)
xpts = range(0,len(doppler))
plt.plot(xpts, doppler)

 
doppler = doppler -signal_freq

plt.figure(2)
xpts = range(0,len(doppler))
plt.plot(xpts, doppler)
plt.show()  
