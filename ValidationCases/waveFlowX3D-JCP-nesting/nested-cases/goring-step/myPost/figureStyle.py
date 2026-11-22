"""作图格式: 取自 waveFlowX-3D-vv/benchmarkSunami/VV01-iswMCC1-flatBottom/myPost/figureStyle.py (serif Times New Roman, 12 pt, 面板标号 (a), (b), ...).
与原文件的区别: text.usetex 只在 LaTeX 真能排版时才开 (本机 PATH 里的 MacPorts latex 缺 type1cm.sty 等宏包, 直接开会整张图失败; 此时改用 /Library/TeX/texbin 的 TeX Live); 关掉时数学符号用 STIX (与 Times 一致), 而不是默认的 DejaVu Sans.
用法: from figureStyle import funApplyStyle, figureLabel; funApplyStyle(); ax.set_title(f'{figureLabel[0]} ...', loc='left')
"""
import os
import random
import shutil

TEXLIVE_BIN = '/Library/TeX/texbin'       # 完整的 TeX Live (macOS); PATH 里的 /opt/local/bin/latex (MacPorts) 缺 type1cm.sty, underscore.sty 等

figureLabel = [f"({chr(ord('a') + i)})" for i in range(26)]


def _renders():
    import matplotlib.pyplot as plt
    old = plt.rcParams['text.usetex']
    try:
        plt.rcParams['text.usetex'] = True
        fig = plt.figure(figsize=(1, 1))
        fig.text(0.5, 0.5, r'probe $x_1$ %d' % random.randrange(10**9))      # 每次不同: matplotlib 按文字缓存 dvi, 相同文字会不调用 latex 就通过
        fig.canvas.draw()
        return True
    except Exception:
        return False
    finally:
        plt.rcParams['text.usetex'] = old
        plt.close('all')


def latexUsable():
    """PATH 里有 latex 不够 (宏包可能不全): 实际排一小段数学文字试一次; 不行再试 TeX Live 的 texbin (成功则把它放到 PATH 最前)."""
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
    """在当前 matplotlib 会话里启用项目统一的格式; 返回 pyplot."""
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
    """文字里的 &: LaTeX 下必须写 \\&, mathtext 下写 & (mathtext 不认 \\&); 在 funApplyStyle() 之后调用."""
    import matplotlib.pyplot as plt
    return r'\&' if plt.rcParams['text.usetex'] else '&'
