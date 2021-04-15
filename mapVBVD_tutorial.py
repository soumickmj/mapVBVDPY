from mapVBVD import mapVBVD
import matplotlib.pyplot as plt
import numpy as np

twix = mapVBVD('meas_MID00051_FID11835_benz_vibe_cap3_qfat_inhale.dat')
# twix = mapVBVD('meas_MID00089_FID13975_sv_ga_5mslice_576sp_fullslice_bo12sp23.dat')

twix = twix[1]

print(twix.image.dataSize[:11])

twix.image.flagDoAverage = True
twix.image.flagRemoveOS = True

print(twix.image.dataSize[:11])

data = twix.image(squeeze=True)

# data = twix.image(shape=[-1, -1, -1, 1,  1,  1,  1,  1,  1,  1, -1], squeeze=False)

fig = plt.figure()
ax = plt.subplot(2, 2, 1)
ax.imshow(np.power(np.abs(data[:, 0, :, 0]), 0.2))

ax = plt.subplot(2, 2, 2)
ax.imshow(np.power(np.abs(data[:, 0, :, 3]), 0.2))

ax = plt.subplot(2, 2, 3)
ax.imshow(np.power(np.abs(data[:, 0, :, 7]), 0.2))

ax = plt.subplot(2, 2, 4)
ax.imshow(np.power(np.abs(data[:, 0, :, 10]), 0.2))

twix.image.flagIgnoreSeg = True
print(twix.image.dataSize[:11])

data = twix.image()
print(data.shape)

plt.show()