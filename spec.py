import tools as tls
import numpy as np
import matplotlib.pyplot as plt
from plot_format import plot_format

#path = "/home/seth/research/wds/wd1235/wd1235-2012-01-25T10:36.fits"
path = "/home/seth/research/wds/wd1235/wd1235-2012-01-30T09:05.fits"
rawwl,rawflux = tls.RawSpectrum(path)

normwl,normflux,normerr = tls.NormNoPlot(path)

plot_format()
plt.plot(rawwl,rawflux)
plt.xlabel("Wavelength [$\AA$]")
plt.ylabel("Flux [$10^{-17}$ erg/s/cm$^{2}$/$\AA$]")
plt.title("KPNO J123549.89+154319.3 Spectrum")
plt.savefig("/home/seth/research/tmp/kpnoSpec.png")

plot_format()
plt.plot(normwl,normflux)
plt.axhline(1.0,color='k',ls='--')
plt.xlabel("Wavelength [$\AA$]")
plt.xlim(3800,5220)
plt.ylabel("Flux [Normalized]")
plt.title("KPNO J123549.89+154319.3 Normalized Spectrum")
plt.savefig("/home/seth/research/tmp/normSpec.png")

