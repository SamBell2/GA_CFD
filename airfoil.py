#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt


class Panel:
    """
    A single panel (flat segment) on the airfoil surface.
    """
    def __init__(self, xa, ya, xb, yb):
        """
        Initialize a panel from point A to point B.

        Parameters:
        -----------
        xa, ya : float
            Coordinates of start point
        xb, yb : float
            Coordinates of end point
        """
        self.xa = xa
        self.ya = ya
        self.xb = xb
        self.yb = yb

        # Calculate center point
        # (control point where we enforce boundary condition)
        self.xc = (xa + xb) / 2
        self.yc = (ya + yb) / 2

        # Calculate panel length
        self.length = np.sqrt((xb - xa)**2 + (yb - ya)**2)

        # Calculate tangent vector (along the panel)
        # Unit vector pointing from A to B
        self.tx = (xb - xa) / self.length
        self.ty = (yb - ya) / self.length

        # Calculate normal vector (perpendicular to panel, pointing outward)
        # For counterclockwise numbering, rotate tangent 90° to the left
        self.nx = -self.ty
        self.ny = self.tx

        # Calculate panel angle (for later use)
        self.beta = np.arctan2(yb - ya, xb - xa)

    def __repr__(self):
        """String representation for printing"""
        return (
            f"Panel: ({self.xa:.3f},{self.ya:.3f}) → \
({self.xb:.3f},{self.yb:.3f}), "
            f"center=({self.xc:.3f},{self.yc:.3f}), length={self.length:.3f}"
        )


class Airfoil:
    def __init__(self,
                 max_camber_dev=0,
                 max_camber_loc=0,
                 thickness=0,
                 num_points=200):
        self.max_camber_dev = max_camber_dev
        self.max_camber_loc = max_camber_loc
        self.thickness = thickness
        self.generateAirfoil(
            max_camber_dev,
            max_camber_loc,
            thickness,
            num_points
        )
        self.name = "NACA ----"

    def __repr__(self):
        self.plot()
        return ''

    def generateAirfoil(self,
                        max_camber_dev,
                        max_camber_loc,
                        thickness,
                        num_points=100) -> None:
        beta = np.linspace(0, np.pi, num_points//2 + 1)
        x = 0.5 * (1 - np.cos(beta))
        yt = (
            5 * thickness
            * (
                0.2969 * np.sqrt(x)
                - 0.1260 * x
                - 0.3516 * x * x
                + 0.2843 * x * x * x
                - 0.1015 * x * x * x * x
            )
        )
        yc = np.zeros_like(x)
        dyc = np.zeros_like(x)
        if max_camber_dev > 0 and max_camber_loc > 0:
            f = x < max_camber_loc
            r = ~f
            yc[f] = max_camber_dev / max_camber_loc**2 * (
                2 * max_camber_loc * x[f] - x[f]**2
            )
            dyc[f] = 2 * max_camber_dev / max_camber_loc**2 * (
                max_camber_loc - x[f]
            )
            yc[r] = max_camber_dev / (1 - max_camber_loc)**2 * (
                (1 - 2 * max_camber_loc) + 2 * max_camber_loc * x[r] - x[r]**2
            )
            dyc[r] = 2 * max_camber_dev / (1 - max_camber_loc)**2 * (
                max_camber_loc - x[r]
            )
        th = np.arctan(dyc)
        xu, yu = x - yt*np.sin(th), yc + yt*np.cos(th)
        xl, yl = x + yt*np.sin(th), yc - yt*np.cos(th)
        X = np.concatenate([xu[::-1], xl[1:]])
        Y = np.concatenate([yu[::-1], yl[1:]])
        self.x_points = X[::-1]
        self.y_points = Y[::-1]

    def plot(self) -> None:
        plt.style.use('seaborn-v0_8-darkgrid')
        plt.switch_backend('TkAgg')
        plt.figure(figsize=(12, 4))
        plt.xlabel("x/c (fraction of chord)", fontsize=12)
        plt.ylabel("y/c (fraction of chord)", fontsize=12)
        plt.legend(fontsize=10)
        plt.grid(True, alpha=0.3)
        plt.axis("equal")  # Keep aspect ratio 1:1 so airfoil isn't distorted
        plt.tight_layout()
        plt.plot(
            self.x_points,
            self.y_points,
            "b-",
            linewidth=2,
            label="Surface"
        )
        plt.title(self.name, fontsize=14, fontweight="bold")
        plt.show()

    def savePlot(self, folder) -> None:
        plt.style.use('seaborn-v0_8-darkgrid')
        plt.switch_backend('TkAgg')
        plt.figure(figsize=(12, 4))
        plt.xlabel("x/c (fraction of chord)", fontsize=12)
        plt.ylabel("y/c (fraction of chord)", fontsize=12)
        plt.legend(fontsize=10)
        plt.grid(True, alpha=0.3)
        plt.axis("equal")  # Keep aspect ratio 1:1 so airfoil isn't distorted
        plt.tight_layout()
        plt.plot(
            self.x_points,
            self.y_points,
            "b-",
            linewidth=2,
            label="Surface"
        )
        plt.title(self.name, fontsize=14, fontweight="bold")
        plt.savefig(folder + '/' + self.name + ".png")

    def savePoints(self, folder) -> None:
        to_save = '\n'.join(
            [f"{x},{y}" for x, y in zip(self.x_points, self.y_points)]
        )
        with open(folder + '/' + self.name + ".csv", 'w') as f:
            f.write(to_save)
