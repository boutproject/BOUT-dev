import numpy as np
import matplotlib.pyplot as plt

its=[0,1,2,3,4,5,6,7,8,9,10,50,100]#,15,20,30,50,100,200]
x = range(12)
y = np.arange(24)/2.0

for it in its:
    b = np.loadtxt("its_"+str(it)+".txt")
    print(b)
    plt.plot(b,'.-',label=it)

b = np.loadtxt("bm.txt")
plt.plot(x,b,'k.:',lw=2)

plt.legend()
plt.savefig("result.pdf")

