"""Oblique linear wave packet in a closed basin with a flat bed: initial fields and the exact linear solution.

Carrier wavelength 1 m (k0 = 2 pi), kh = 1 (h = 0.159155 m), amplitude a = 1e-3 h, direction THETA from the x axis,
isotropic Gaussian envelope of standard deviation SIGMA centred at X0 at t = 0. The packet is the real part of the
complex field a exp(-|x - x0|^2 / 2 sigma^2) exp(i k0 . (x - x0)); every Fourier component of that field travels in the
direction of its own wave vector with the linear dispersion relation (the components near -k0 that a real cosine
would add are below 1e-30 of the peak for k0 sigma = 2 pi). u, v and w of each component are those of linear theory
at the sigma-layer centres of the still-water column.
exact_eta evaluates the same solution at any (x, y, t) from the components above 1e-13 of the peak; it describes an
unbounded basin, so the walls must stay unreached (checked by wall_level).
"""
import math
import numpy as np

GRAV, LAMBDA = 9.81, 1.0
K0 = 2 * math.pi / LAMBDA
KH = 1.0
H = KH / K0
AMP = 1e-3 * H
THETA = math.radians(30.0)
SIGMA = 1.0
LX, LY = 25.0, 21.0
X0 = np.array([6.0, 6.0])
TRAVEL = 10.0                                   # path length of the envelope centre at T_END
CHILD = (8.3, 12.3, 6.5, 10.5)                  # x0, x1, y0, y1: 4 x 4 m centred on the path at s = 5 m
DIR = np.array([math.cos(THETA), math.sin(THETA)])


def omega(k):
    return np.sqrt(GRAV * k * np.tanh(k * H))


def group_velocity():
    return 0.5 * float(omega(K0)) / K0 * (1 + 2 * KH / math.sinh(2 * KH))


T_END = TRAVEL / group_velocity()


def path(s):
    return X0 + s * DIR


def complex_packet(x, y):
    dx, dy = x - X0[0], y - X0[1]
    return AMP * np.exp(-(dx ** 2 + dy ** 2) / (2 * SIGMA ** 2)) * np.exp(1j * K0 * (DIR[0] * dx + DIR[1] * dy))


def _grid(nx, ny):
    dx, dy = LX / nx, LY / ny
    return (np.arange(nx) + .5) * dx, (np.arange(ny) + .5) * dy, dx, dy


def initial_fields(nx, ny, kz):
    """eta (ny, nx); u, v, w (kz, ny, nx) at the sigma-layer centres of the still-water column."""
    x, y, dx, dy = _grid(nx, ny)
    X, Y = np.meshgrid(x, y)
    spec = np.fft.fft2(complex_packet(X, Y))
    kx = 2 * math.pi * np.fft.fftfreq(nx, dx)
    ky = 2 * math.pi * np.fft.fftfreq(ny, dy)
    KX, KY = np.meshgrid(kx, ky)
    K = np.hypot(KX, KY)
    K[0, 0] = 1.0
    spec[0, 0] = 0.0
    W = omega(K)
    z = -H + (np.arange(kz) + .5) / kz * H
    den = 1 - np.exp(-2 * K * H)
    u, v, w = (np.empty((kz, ny, nx)) for _ in range(3))
    for j, zj in enumerate(z):
        ch = (np.exp(K * zj) + np.exp(-K * (zj + 2 * H))) / den      # cosh(K (z + h)) / sinh(K h)
        sh = (np.exp(K * zj) - np.exp(-K * (zj + 2 * H))) / den      # sinh(K (z + h)) / sinh(K h)
        u[j] = np.real(np.fft.ifft2(spec * W * ch * KX / K))
        v[j] = np.real(np.fft.ifft2(spec * W * ch * KY / K))
        w[j] = np.real(np.fft.ifft2(-1j * spec * W * sh))
    return np.real(np.fft.ifft2(spec)), u, v, w


class Exact:
    """Exact linear eta(x, y, t) from the significant components of a fine spectral representation."""

    def __init__(self, nfine=(1000, 840), tol=1e-13):
        nx, ny = nfine
        x, y, dx, dy = _grid(nx, ny)
        X, Y = np.meshgrid(x, y)
        spec = np.fft.fft2(complex_packet(X, Y)) / (nx * ny)
        kx = 2 * math.pi * np.fft.fftfreq(nx, dx)
        ky = 2 * math.pi * np.fft.fftfreq(ny, dy)
        KX, KY = np.meshgrid(kx, ky)
        keep = np.abs(spec) > tol * np.abs(spec).max()
        self.c, self.kx, self.ky = spec[keep], KX[keep], KY[keep]
        self.w = omega(np.hypot(self.kx, self.ky))
        self.x0, self.y0 = x[0], y[0]

    def eta(self, x, y, t):
        """x, y scalars or arrays of one shape; t an array: returns shape t.shape + x.shape."""
        x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
        t = np.asarray(t, float)
        sp = np.exp(1j * (np.multiply.outer(x - self.x0, self.kx) + np.multiply.outer(y - self.y0, self.ky)))
        out = np.empty(t.shape + x.shape)
        for i, ti in np.ndenumerate(t):
            out[i] = np.real(sp @ (self.c * np.exp(-1j * self.w * ti)))
        return out


def wall_level(nx=500, ny=420, times=None):
    """Largest |eta| / a on the four walls over the run (unbounded solution), for the margin check."""
    ex = Exact()
    x, y, _, _ = _grid(nx, ny)
    xs = np.r_[x, x, np.zeros(ny), np.full(ny, LX)]
    ys = np.r_[np.zeros(nx), np.full(nx, LY), y, y]
    times = np.linspace(0, T_END, 25) if times is None else times
    return float(np.abs(ex.eta(xs, ys, times)).max() / AMP)
