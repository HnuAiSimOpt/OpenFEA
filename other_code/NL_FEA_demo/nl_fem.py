# -------------------------------------------------------------------------
# Copyright:  WH team
# Author: YinJichao <jichaoyinyjc@163.com>
# Completion date:  XXX
# Description: XXX
# -------------------------------------------------------------------------
from read_inp import *
from fem_geo_nl import *


path = ".\\geo_model\\2d\\Job-NL.inp"

if __name__ == "__main__":
    # read
    node, element, load, fixed, em, nu, load_axis, load_value = read_inp(path)
    node = np.delete(node, 0, axis=1)
    element = np.delete(element, 0, axis=1)
    # analysis
    nl_fem = GeoNL3D()
    nl_fem.geo_nonlinear_analysis(node, element, load, fixed, em, nu, load_axis, load_value)


