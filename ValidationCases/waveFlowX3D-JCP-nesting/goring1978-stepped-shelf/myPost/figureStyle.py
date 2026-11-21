"""Figure style: serif Times New Roman, 12 pt, panel labels (a), (b), ...
LaTeX text rendering is switched on only if a working LaTeX is found; otherwise mathtext with STIX fonts is used.
Usage: from figureStyle import funApplyStyle, figureLabel; plt = funApplyStyle(); ax.set_title(f'{figureLabel[0]} ...', loc='left')
"""
import os
import random
import shutil

TEXLIVE_BIN = '/Library/TeX/texbin'       # full TeX Live on macOS, tried when the latex on PATH lacks packages

figureLabel = [f"({chr(ord('a') + i)})" for i in range(26)]


def _renders():
    import matplotlib.pyplot as plt
    old = plt.rcParams['text.usetex']
    try:
        plt.rcParams['text.usetex'] = True
        fig = plt.figure(figsize=(1, 1))
        fig.text(0.5, 0.5, r'probe $x_1$ %d' % random.randrange(10**9))      # a new text each time: matplotlib caches dvi files by text
        fig.canvas.draw()
        return True
    except Exception:
        return False
    finally:
        plt.rcParams['text.usetex'] = old
        plt.close('all')


def latexUsable():
    """A latex on PATH is not enough (packages may be missing): typeset a short math text once; if that fails, try
    the TeX Live texbin (and put it first on PATH if it works)."""
    if shutil.which('latex') is not None and _renders():
        return True
    if os.path.isdir(TEXLIVE_BIN) and os.path.exists(os.path.join(TEXLIVE_BIN, 'latex')):
        old = os.environ['PATH']
        os.environ['PATH'] = TEXLIVE_BIN + os.pathsep + old
        if _renders():
            return True
        os.environ['PATH'] = old
    return False


def funApplyStyle():
    """Apply the figure style to the current matplotlib session; return pyplot."""
    import matplotlib.pyplot as plt
    usetex = latexUsable()
    plt.rcParams.update({
        'font.family': 'serif',
        'font.serif': ['Times New Roman', 'STIXGeneral', 'DejaVu Serif'],
        'font.size': 12.0,
        'text.usetex': usetex,
        'mathtext.fontset': 'stix',
        'axes.unicode_minus': False,
    })
    return plt


def amp():
    """An ampersand in text: \\& under LaTeX, & under mathtext (which does not accept \\&); call after funApplyStyle()."""
    import matplotlib.pyplot as plt
    return r'\&' if plt.rcParams['text.usetex'] else '&'
