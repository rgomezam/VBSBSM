# Import Library

import matplotlib.pyplot as plt

# Define data coordinates

x = [0, 1, 2, 3, 4, 5]
y = [1.5, 3, 5.3, 6, 10, 2]
y2 = [3, 5,  6,  5,  9, 13]

print('y2= ', y2)

# Plot

plt.plot( y, y2, '-o', color='red')

# fill_between

plt.fill_between(x, y)

# Add title

plt.suptitle('Simple fill_between Example', fontweight='bold')

# Display

plt.show()