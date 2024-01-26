# -------------------------------------------------------------------------
# Copyright:  WH team
# Author: YinJichao <jichaoyinyjc@163.com>
# Completion date:  XXX
# Description: XXX
# -------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def show_dis(x, y, u_x, u_y, save_path):
    X = np.expand_dims(x + u_x, axis=1)
    Y = np.expand_dims(y + u_y, axis=1)
    plt.rcParams['font.family'] = 'Times New Roman'
    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(16., 6.))
    # subplot 1
    cfig1 = axs[0].scatter(X, Y, c=u_x, vmin=u_x.min(), vmax=u_x.max(), cmap='coolwarm')
    axes_1 = inset_axes(
        axs[0],
        width=0.15,  # width: 5% of parent_bbox width
        height="100%",  # height: 50%
        loc="lower left",
        bbox_to_anchor=(1.05, 0., 1, 1),
        bbox_transform=axs[0].transAxes,
        borderpad=0,
    )
    axs[0].set_aspect('equal')
    axs[0].set_title('displacement X')
    plt.colorbar(cfig1, cax=axes_1)
    # subplot 2
    cfig2 = axs[1].scatter(X, Y, c=u_y, marker='o', vmin=u_y.min(), vmax=u_y.max(), cmap='coolwarm')
    axes_2 = inset_axes(
        axs[1],
        width=0.15,  # width: 5% of parent_bbox width
        height="100%",  # height: 50%
        loc="lower left",
        bbox_to_anchor=(1.05, 0., 1, 1),
        bbox_transform=axs[1].transAxes,
        borderpad=0,
    )
    axs[1].set_aspect('equal')
    axs[1].set_title('displacement Y')
    plt.colorbar(cfig2, cax=axes_2)
    # save
    plt.savefig(save_path, dpi=600, bbox_inches='tight')
    plt.close()
