import numpy as np
import matplotlib.pyplot as plt 

#h = np.array([])


#h=np.append(h, []) 

h0 = 1.000000000000000000E-01    
h1 = 1.000000000000000000E-02
h2 = 1.000000000000000000E-03    
h3 = 1.000000000000000000E-04    
h4 = 1.000000000000000000E-05    
h5 = 1.000000000000000000E-06    
h6 = 1.000000000000000000E-07    
h = np.array([h0,h1,h2,h3,h4,h5,h6])

en0= 1.110223024625157000E-16
en1= 1.998401444325282000E-15
en2= 1.226796442210798000E-13
en3= 5.639155808978558000E-12
en4= 6.027178756085050000E-11
en5= 2.634219269381788000E-07
en6= 1.602743389006456000E-06    
En = np.array([en0,en1,en2,en3,en4,en5,en6])

plt.subplot(221)
plt.plot(h,En)
plt.yscale('linear')
plt.title('En x h')
plt.grid(True)

plt.subplot(222)
plt.plot(h,En)
plt.yscale('log')
plt.xscale('log')
plt.title('En x h')
plt.grid(True)


plt.show()
