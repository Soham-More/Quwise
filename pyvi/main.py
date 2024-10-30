from pyvisual import *

pyvi = PyVi('v.pyvi')

#pyvi.plot_iteration('potential', 1, label_axes=True, plt_format='')

#pyvi.plot_common('position', 1)

i = 0

pyvi.display_all_sections()

#pyvi.plot_sections('position', ['potential', 'Fermi lvl'], i)

#pyvi.plot_sections('position', ['charge', 'charge derivative'], i)

#pyvi.plot_sections('position', ['residual', 'Jacobian'], i)

#pyvi.plot_iteration('residual', i, label_axes=True, plt_format='')

plt.show()

