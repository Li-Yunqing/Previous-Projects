import numpy as np 
import matplotlib.pyplot as plt 

a_values = np.arange(2, 4, 0.001)
steps = 1000
initial_value = 0.6

x_values = np.zeros((steps, a_values.__len__()))
x_values[0] = initial_value

#循环这段很不好,没有并行化，要重写
i = 0
for a in a_values:
    for step in range(0, steps-1):
        x_values[step+1, i] = a * x_values[step, i] * (1 - x_values[step, i])
    i = i + 1    

#plot这段也要改，自由度还不够高
plt.figure(1)
i = 0
for x in x_values[100:,:]:
    plt.scatter(a_values, x)
    i = i + 1



initial_value = 0.60000001
a = 3.947
b = np.zeros(steps)
b[0] = initial_value
for step in range(0, steps-1):
    b[step+1] = a * b[step] * (1 - b[step])

plt.figure(2)
plt.plot(range(0,50), x_values[:50, 1947])
plt.plot(range(0,50), b[:50])





plt.figure(2)
plt.subplot(231)
plt.plot(range(950,1000), x_values[950:, 500])

plt.subplot(232)
plt.plot(range(950,1000), x_values[950:, 1250])

plt.subplot(233)
plt.plot(range(950,1000), x_values[950:, 1500])

plt.subplot(234)
plt.plot(range(950,1000), x_values[950:, 1560])

plt.subplot(235)
plt.plot(range(950,1000), x_values[950:, 1947])

plt.subplot(236)
plt.plot(range(950,1000), x_values[950:, 1847])

plt.show()

