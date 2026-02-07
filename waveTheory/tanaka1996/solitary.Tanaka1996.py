import matplotlib

# 根据系统环境自动选择后端，防止在无显示器服务器上报错
try:
    matplotlib.use('TkAgg')
except:
    matplotlib.use('Agg')

import numpy as np
from scipy.optimize import fsolve, brentq
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import griddata, interp1d  # 新增插值工具
import matplotlib.pyplot as plt
import os

# 设置绘图字体
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman"],
    "font.size": 14.0,
    "text.usetex": True  # 如有LaTeX环境可设为True
})

# 环境变量设置，修复部分环境下的distutils报错
os.environ['SETUPTOOLS_USE_DISTUTILS'] = 'stdlib'


class TanakaSolver:
    """
    Tanaka (1986) 孤立波全非线性求解器 (v8: 物理域调整版 [-1, eta])

    ------------------------------------------------------------------
    【物理模型与数学原理】

    1. 控制方程:
       考虑二维、无粘、不可压缩、无旋的孤立水波。
       引入复速度势 (Complex Velocity Potential) W = Phi + i*Psi。
       其中:
         - Phi (Velocity Potential): 速度势函数
         - Psi (Stream Function): 流函数

       物理域被共形映射到复平面上的一个无限长条带 (Strip Domain):
         - 0 < Psi < 1
         - Psi = 1: 自由表面 (Free Surface)
         - Psi = 0: 平坦底部 (Flat Bottom)

    2. 核心变量 (Log-Hodograph Variables):
       引入对数复速度函数 Omega = ln(dW/dz) = tau - i*theta
       其中:
         - q = |dW/dz|: 速度大小
         - tau = ln(q): 对数速度
         - theta: 速度矢量与水平轴的夹角 (流线倾角)
         - z = x + i*y: 物理坐标

    3. 边界条件:
       (1) 底部条件 (Psi=0): theta = 0 (不可穿透)
       (2) 自由表面动力学条件 (Psi=1): 伯努利方程 (Bernoulli Equation)
           q^3 = qc^3 - (3/F^2) * int_{0}^{Phi} sin(theta) dPhi  [Tanaka Eq. 4]
           其中 F 是弗劳德数 (Froude Number), qc 是波峰速度。

    4. 求解策略 (Tanaka Algorithm):
       利用边界积分方程 (Boundary Integral Equation) 将域内问题转化为边界(自由表面)上的非线性积分方程组，
       并通过 Picard 迭代求解。
    ------------------------------------------------------------------
    """

    def __init__(self, eps=None, alpha=0.0125, m=5, qc=0.11, tol=1e-7, phi_max=40, N=30):
        """
        初始化求解器参数。

        参数详解:
        ----------------
        eps : float (Optional)
            目标相对波高 H/h (Wave Height / Water Depth).
            若提供此参数，程序将自动通过迭代搜索找到对应的输入参数 qc。

        alpha, m : float, int
            坐标变换参数 [Tanaka Eq. 6]: Phi = alpha * gamma + gamma^m
            物理意义:
              孤立波在波峰处 (Phi=0) 变化最剧烈。该变换将计算网格 gamma 映射到物理势 Phi。
              在 gamma 均匀分布的情况下，Phi 在 0 附近会被显著加密，从而提高波峰处的计算精度。
              - alpha 控制波峰附近的网格密度。
              - m 控制远场的拉伸程度。

        qc : float
            波峰处的流体速度 (Wave Crest Velocity). 无量纲化单位: c (波速).
            物理意义:
              在随波坐标系 (Wave Frame) 下，波形静止，水流向后流动。
              波峰处水深变浅，流速减小。
              - qc = q_crest / c
              - qc -> 1: 对应极小振幅波 (线性波).
              - qc -> 0: 对应极限波高 (120度尖角, 流速停滞点).

        phi_max : float
            计算域在势空间 Phi 上的截断范围 [-phi_max, phi_max].
            物理意义: 模拟无限远边界条件 (Phi -> infinity)。

        N : int
            计算网格点数 (正半轴)。总节点数为 2N-1。

        tol : float
            迭代收敛容差 (针对 F^2)。
        """
        # 1. 初始化配置参数
        self.alpha = alpha
        self.m = m
        self.tol = tol
        self.phi_max = phi_max
        self.N = N
        self.interp_order = 8  # 使用8点拉格朗日插值

        # 2. 声明所有全局变量 (Initialize all global variables to None)
        self.qc = qc
        self.h = None
        self.H = None
        self.L0 = None
        self.L00 = None

        # 计算结果存储变量
        self.F2 = None  # 弗劳德数平方 (无量纲波速平方)
        self.Mass = None  # 质量
        self.PotentialEn = None  # 势能
        self.Impulse = None  # 冲量
        self.KineticEn = None  # 动能

        # 网格与场变量 (中间变量)
        self.r_full = None  # Gamma 网格 (计算域)
        self.phi_full = None  # 势函数 Phi 网格
        self.theta_full = None  # 表面流线倾角
        self.last_tau_right = None  # 表面对数速度

        # 物理空间数据 (表面)
        self.grid_x = None
        self.data_eta = None
        self.data_us = None
        self.data_ws = None

        # 物理空间数据 (内部/中心线)
        self.grid_z = None
        self.data_uc = None
        self.data_wc = None

        # 物理空间数据 (Sigma网格结果) [新增]
        self.data_U = None
        self.data_W = None
        self.grid_X = None  # 对应的物理坐标 X 网格
        self.grid_Z = None  # 对应的物理坐标 Y 网格
        self.grid_Sigma = None

        # 积分权重矩阵 (加速计算)
        self.W_int_std = None
        self.W_interp_std = None

        # 3. 预计算权重
        self._precompute_weights()

        # 4. 如果提供了 eps，执行完整计算流程并赋值给 self 变量
        if eps is not None:
            # 反算 qc (利用 Brent 方法寻找对应的波峰速度)
            print(f"[初始化] 目标波高 H/h = {eps}，正在反算对应的 qc 参数...")
            self.qc = self._find_qc_from_height(eps)
            print(f"[初始化] 反算完成: qc = {self.qc:.8f}")

            self.h = 1.0
            self.H = eps
            self.eps = self.H / self.h
            self.L0 = 3.5 * self.h / np.sqrt(self.H / self.h)
            self.L00 = np.ceil(self.L0)

        if True:
            # (1) 求解流场 (Pure Function Call)
            # 核心算法入口
            r_vec, theta_pos, F2, tau_right, theta_full, r_full, phi_full = self.solve(self.qc)

            # 存储流场核心结果
            # self.r_vec = r_vec
            # self.theta_pos = theta_pos
            self.F2 = F2
            self.last_tau_right = tau_right
            self.theta_full = theta_full
            self.r_full = r_full
            self.phi_full = phi_full

            # ... 在 TanakaSolver.__init__ 调用计算部分 ...
            if True:
                # (2) 计算物理变量
                grid_x, eta, us, ws, M, C, V, I, K = self.compute_physical_variables(
                    r_vec, theta_pos, F2, tau_right
                )

                # 存储结果
                self.grid_x = grid_x
                self.data_eta = eta
                self.data_us = us
                self.data_ws = ws
                self.Mass = M
                self.Circulation = C  # 存储环量 [cite: 312, 1294]
                self.PotentialEn = V
                self.Impulse = I
                self.KineticEn = K

        if True:
            # (3) 计算内部速度 (Pure Function Call)
            # 利用解析延拓计算波峰下方的速度分布
            psi_vec = np.linspace(1, 0, 50)
            x_shift, grid_z, uc, wc = self.compute_internal_velocity_at_phi(
                0.0, psi_vec, theta_full, phi_full, r_full
            )

            # 存储内部速度结果
            self.grid_z = grid_z
            self.data_uc = uc
            self.data_wc = wc

        if True:
            # --- 步骤 1: 利用对称性构建全域表面数据 ---
            # 表面 x 坐标是奇函数: x(-phi) = -x(phi)
            full_x = np.concatenate((-self.grid_x[::-1][:-1], self.grid_x))
            # 表面位移 eta 是偶函数: eta(-phi) = eta(phi)
            full_eta = np.concatenate((self.data_eta[::-1][:-1], self.data_eta))
            # 水平速度 us 是偶函数: us(-phi) = us(phi)
            full_us = np.concatenate((self.data_us[::-1][:-1], self.data_us))
            # 垂直速度 ws 是奇函数: ws(-phi) = -ws(phi)
            full_ws = np.concatenate((-self.data_ws[::-1][:-1], self.data_ws))

            # 更新类成员变量为全域数据
            self.grid_x = full_x
            self.data_eta = full_eta
            self.data_us = full_us
            self.data_ws = full_ws

            # --- 步骤 2: 初始化全网格 ---
            self.grid_phi_full = self.phi_full
            self.grid_psi = psi_vec

            # 初始化矩阵
            N_phi = len(self.phi_full)
            N_psi = len(psi_vec)
            self.grid_X = np.zeros((N_phi, N_psi))
            self.grid_Z = np.zeros((N_phi, N_psi))
            self.grid_Sigma = np.zeros((N_phi, N_psi))
            self.data_U = np.zeros((N_phi, N_psi))
            self.data_W = np.zeros((N_phi, N_psi))

            # 建立 Phi 到 Surface X 的映射
            f_phi_to_x_surf = interp1d(self.phi_full, full_x, bounds_error=False, fill_value="extrapolate")

            # --- 步骤 3: 单次遍历计算全流场 ---
            for index, phi in enumerate(self.phi_full):
                # 计算剖面
                x_shift, grid_z, tmp_u, tmp_w = self.compute_internal_velocity_at_phi(
                    phi, psi_vec, theta_full, phi_full, r_full
                )

                # 强制边界一致性 (Index 0 is surface)
                tmp_u[0] = full_us[index]
                tmp_w[0] = full_ws[index]
                grid_z[0] = full_eta[index]

                # 计算 Sigma
                tmp_sigma = (grid_z + 1.0) / (full_eta[index] + 1.0)

                # 计算绝对 X 坐标 (表面 X + 内部偏移)
                surface_x = f_phi_to_x_surf(phi)
                real_x = surface_x + x_shift

                # 填充矩阵
                self.grid_X[index, :] = real_x
                self.grid_Z[index, :] = grid_z
                self.grid_Sigma[index, :] = tmp_sigma
                self.data_U[index, :] = tmp_u
                self.data_W[index, :] = tmp_w

        



    def _solve_r(self):
        """辅助函数: 求解计算域边界 r_max (对应 Phi_max)"""
        # 方程: alpha * r + r^m = phi_max
        func = lambda r: self.alpha * r + r ** self.m - self.phi_max
        r0 = (self.phi_max / self.alpha) ** (1 / self.m)
        r_max = fsolve(func, r0)[0]
        return r_max

    def _precompute_weights(self):
        """
        预计算权重矩阵。
        原理: 将高阶拉格朗日积分和插值转化为矩阵乘法，将 O(N^2) 复杂度降为 O(N)。
        """
        n = self.interp_order
        std_nodes = np.arange(n)
        self.W_int_std = np.zeros((n - 1, n))
        self.W_interp_std = np.zeros((n, n))
        for j in range(n):
            y_basis = np.zeros(n)
            y_basis[j] = 1.0
            poly = np.polynomial.Polynomial.fit(std_nodes, y_basis, n - 1)
            poly_int = poly.integ()

            # 积分权重
            for k in range(n - 1):
                self.W_int_std[k, j] = poly_int(k + 1) - poly_int(k)

            # 中点插值权重
            for k in range(n):
                self.W_interp_std[k, j] = poly(k + 0.5)

    def _lagrange_int_all_optimized(self, nodes, f):
        """数值积分函数: 基于预计算权重的快速积分"""
        n = self.interp_order
        N_total = len(nodes)
        num_sub = (N_total - 1) // (n - 1)
        integration_set = np.zeros(N_total)
        dr = nodes[1] - nodes[0]
        real_W = self.W_int_std * dr

        for k in range(num_sub):
            idx_start = k * (n - 1)
            idx_end = (k + 1) * (n - 1) + 1
            integration_set[idx_start + 1:idx_end] = real_W @ f[idx_start:idx_end]

        st = num_sub * (n - 1)
        en = N_total
        if st < en - 1:
            rel_st = st - (en - n)
            integration_set[st + 1:en] = real_W[rel_st:, :] @ f[en - n:en]

        return np.cumsum(integration_set)

    def _midpoint_calculate_optimized(self, nodes, f):
        """插值函数: 将节点值插值到网格中点"""
        n = self.interp_order
        N_total = len(nodes)
        num_sub = (N_total - 1) // (n - 1)
        mid_f = np.zeros(N_total)
        for k in range(num_sub):
            idx_start = k * (n - 1)
            idx_end = (k + 1) * (n - 1) + 1
            mid_f[idx_start:idx_end - 1] = self.W_interp_std[:n - 1, :] @ f[idx_start:idx_end]
        st = num_sub * (n - 1)
        en = N_total
        if st < en:
            mid_f[st:en] = self.W_interp_std[st - (en - n):en - (en - n), :] @ f[en - n:en]
        return mid_f

    def _compute_kernel_matrix(self, xvec, yvec):
        """
        计算边界积分方程的核矩阵 (Kernel Matrix)。

        数学原理:
        theta(Phi) = - integral { K(Phi, Phi') * tau(Phi') } dPhi'
        核函数 K(Phi, Phi') = 1 / [ 2 * sinh( pi/2 * (Phi' - Phi) ) ]

        技术细节:
        采用交错网格 (xvec 是节点, yvec 是中点) 避免处理 Phi' = Phi 处的奇异性。
        """
        X = xvec[:, None]
        Y = yvec[None, :]

        # 坐标变换 gamma -> Phi
        phi_x = self.alpha * X + X ** self.m
        phi_y = self.alpha * Y + Y ** self.m

        Mat = 1.0 / (2.0 * np.sinh((phi_y - phi_x) * np.pi / 2.0))
        # 雅可比行列式 dPhi/dgamma
        jacobian = self.alpha + self.m * Y ** (self.m - 1)

        # 梯形积分权重
        W_trap = np.zeros_like(yvec)
        if len(yvec) > 1:
            W_trap[0] = (yvec[1] - yvec[0]) / 2.0
            W_trap[-1] = (yvec[-1] - yvec[-2]) / 2.0
            W_trap[1:-1] = (yvec[2:] - yvec[:-2]) / 2.0

        return Mat * jacobian * W_trap

    def solve(self, qc_val):
        """
        主求解过程：执行 Picard 迭代。

        算法流程:
        1. 猜测初始 tau (对数速度)。
        2. 利用边界积分方程计算 theta (流向角)。
        3. 利用伯努利方程更新 F^2 (波速) 和 q (速度大小)。
        4. 迭代直到 F^2 收敛。

        Returns values ONLY. Does not modify self.
        """
        # --- 1. 网格生成 ---
        r_max = self._solve_r()
        dr = r_max / self.N
        r_vec = dr * np.arange(self.N)

        # 全域对称网格 (用于积分)
        r = np.concatenate((-r_vec[1:][::-1], r_vec))
        # Equantion (6)
        phi_vec_posi = self.alpha * r_vec + r_vec ** self.m

        # 中点网格 (用于源点积分，避开奇异性)
        r_vec_half = dr * (np.arange(self.N) + 0.5)
        midr = np.concatenate((-r_vec_half[::-1], r_vec_half))

        # --- 2. 初始猜测 ---
        # 经验公式: 假设速度呈指数衰减
        # tau = ln(q), tau(0) = ln(q_c)
        tau_func = lambda x: np.log(qc_val) * np.exp(-0.5 * np.abs(x))
        tau_right = tau_func(phi_vec_posi)
        tau_left = tau_func(-phi_vec_posi[::-1])

        F2 = 0.0
        F2_old = 1e15
        iteration = 0

        # 声明 theta 变量，防止未收敛时报错
        theta = np.zeros_like(r)
        theta_pos = np.zeros_like(r_vec)

        # --- 3. 迭代循环 (Picard Iteration) ---
        while abs(F2_old - F2) > self.tol and iteration < 100:
            iteration += 1

            # A. 准备积分源项: 将 tau 插值到中点
            trans_tau_left = tau_left[::-1]
            mid_tau = np.concatenate((
                self._midpoint_calculate_optimized(r_vec, trans_tau_left)[::-1],
                self._midpoint_calculate_optimized(r_vec, tau_right)
            ))

            # B. 求解边界积分方程 [Tanaka Eq. 3]
            # theta = - Kernel * tau
            M = self._compute_kernel_matrix(r, midr)
            theta = - M @ mid_tau
            theta_pos = theta[self.N - 1:]

            # C. 更新波速参数 F^2 [Tanaka Eq. 5]
            # 物理条件: 无穷远处流速为 1 (Solvability Condition)
            # F^2 = 3 * Int(sin(theta)) / (qc^3 - 1)
            dr_vals = r_vec[1:] - r_vec[:-1]
            w_vec = np.zeros(len(r_vec))
            w_vec[0] = dr_vals[0] / 2.0
            w_vec[-1] = dr_vals[-1] / 2.0
            w_vec[1:-1] = (dr_vals[:-1] + dr_vals[1:]) / 2.0

            integrand_F2 = np.sin(theta_pos) * (self.alpha + self.m * r_vec ** (self.m - 1))
            integral_val = np.dot(integrand_F2, w_vec)

            F2 = 3 * integral_val / (qc_val ** 3 - 1)

            # D. 更新流场速度 q 和 tau [Tanaka Eq. 4]
            # 物理条件: 自由表面伯努利方程 (Bernoulli Equation)
            # q^3 = qc^3 - (3/F^2) * Int(sin(theta))
            integral_vec_right = self._lagrange_int_all_optimized(r_vec, integrand_F2)

            theta_neg_reversed = theta[self.N - 1::-1]
            integrand_neg = np.sin(theta_neg_reversed) * (self.alpha + self.m * r_vec ** (self.m - 1))
            trans_integral_vec_left = -self._lagrange_int_all_optimized(r_vec, integrand_neg)

            # 更新 q^3，并限制最小值防止数值错误
            new_q3_right = np.maximum(qc_val ** 3 - (3.0 / F2) * integral_vec_right, 1e-10)
            new_q3_left = np.maximum(qc_val ** 3 - (3.0 / F2) * trans_integral_vec_left, 1e-10)

            # 更新 tau = ln(q)
            tau_right = np.log(new_q3_right) / 3.0
            tau_left = (np.log(new_q3_left) / 3.0)[::-1]

        # 构造辅助全域变量返回
        phi_full = self.alpha * r + r ** self.m if (self.m % 2 != 0) else self.alpha * r + np.sign(r) * np.abs(
            r) ** self.m
        r_full = r
        theta_full = theta

        return r_vec, theta_pos, F2, tau_right, theta_full, r_full, phi_full

    def _find_qc_from_height(self, eps):
        """
        私有方法: 反算 qc。
        给定目标波高 eps (H/h)，寻找对应的输入参数 qc (波峰速度)。
        """

        def objective(qc_guess):
            if not (0 < qc_guess < 1): return -1e9
            # 实例化一个临时求解器，不带 eps 防止递归
            solver = self.__class__(eps=None, qc=qc_guess, N=30, tol=1e-5)

            # 显式调用纯计算方法，获取波高
            r, theta, F2, tau, _, _, _ = solver.solve(qc_guess)
            _, eta, _, _, _, _, _, _, _  = solver.compute_physical_variables(r, theta, F2, tau)

            # 误差 = 计算波高 - 目标波高
            return (eta[0] - 0.0) - eps

        if eps > 0.833:
            print("警告: 目标波高 H/h > 0.833，超出孤立波物理极限。")
            return 0.01

        return brentq(objective, 0.001, 0.99)

    def compute_physical_variables(self, r_vec, theta_pos, F2, tau_right):
        """
        计算波形、表面速度及积分物理量 (质量、环量、动能、势能等)。
        """
        tau = tau_right
        q = np.exp(tau)
        c = np.sqrt(F2)
        # dphi_dr, Eq(6)
        jacobian = self.alpha + self.m * r_vec ** (self.m - 1)

        # 1. 计算坐标微分项
        dx_dphi = (np.cos(theta_pos) / q)
        dz_dphi = (np.sin(theta_pos) / q)

        # dx, dz 是关于 r 的微分
        dx = dx_dphi * jacobian
        dz = dz_dphi * jacobian

        # 2. 积分得到物理坐标
        x = cumulative_trapezoid(dx, r_vec, initial=0)
        z = cumulative_trapezoid(dz, r_vec, initial=0)
        z_abs = z - z[-1]  # 调整使无穷远处 y=0
        eta = z_abs

        # 3. 计算积分物理量 (利用对称性计算半域再乘以2)
        # 质量 M = 2 * \int eta dx
        M_half = cumulative_trapezoid(eta, x)[-1]

        # 势能 V = 2 * \int 0.5 * eta^2 dx
        V_half = 0.5 * cumulative_trapezoid(eta ** 2, x)[-1]

        # 环量 C = 2 * \int (cos(theta)/q - 1) dPhi
        # 注意: 这里的 q 是无量纲速度 q/c
        integrand_C = (np.cos(theta_pos) / q - 1.0) * jacobian
        C_half = cumulative_trapezoid(integrand_C, r_vec, initial=0)[-1]

        Mass = 2 * M_half
        PotentialEn = 2 * V_half
        Circulation = 2 * C_half

        Impulse = c * Mass
        # 利用恒等式计算动能 K
        KineticEn = 0.5 * (c * Impulse - 2 * PotentialEn)

        # 返回所有计算值
        u_wave = q * np.cos(theta_pos)
        w_wave = q * np.sin(theta_pos)
        return x, z_abs, (1 - u_wave)*c, -w_wave * c, Mass, Circulation, PotentialEn, Impulse, KineticEn

    def compute_internal_velocity_at_phi(self, target_phi, psi_vec, theta_full, phi_full, r_full):
        """
        计算内部流场速度及物理坐标 (修正版：修复水平位移符号错误)。
        """
        dr = r_full[1] - r_full[0]
        jacobian = self.alpha + self.m * np.abs(r_full) ** (self.m - 1)

        u_list, v_list = [], []

        # 1. 解析延拓计算速度
        for psi in psi_vec:
            W = target_phi + 1j * psi
            arg = np.pi * (W - phi_full)
            kernel = 0.5 * (1.0 - np.tanh(arg / 2.0))

            integrand = theta_full * kernel * jacobian
            omega_val = np.sum(integrand) * dr

            tau_val = omega_val.real
            theta_val = -omega_val.imag
            q = np.exp(tau_val)

            # u_fixed, v_fixed 是实验室坐标系下的速度分量
            # v_fixed = q * sin(theta)
            u_fixed = 1.0 - q * np.cos(theta_val)
            v_fixed = q * np.sin(theta_val)

            u_list.append(u_fixed)
            v_list.append(v_fixed)

        u_arr = np.array(u_list)
        v_arr = np.array(v_list)

        # 2. 积分物理坐标
        # q_wave^2 = q^2
        q_wave_sq = (1.0 - u_arr) ** 2 + v_arr ** 2

        # dy = (cos(theta)/q) * dPsi
        dy_dpsi = (1.0 - u_arr) / q_wave_sq

        # dx = -(sin(theta)/q) * dPsi
        # v_arr = q * sin(theta)
        # 所以 -sin(theta)/q = - (v_arr/q) / q = -v_arr / q^2
        dx_dpsi = -v_arr / q_wave_sq  # <--- 【关键修正】这里必须加负号

        # 从表面 (psi=1) 向下积分 (psi_vec 递减，dPsi 为负)
        y_from_surf = cumulative_trapezoid(dy_dpsi, psi_vec, initial=0)
        x_from_surf = cumulative_trapezoid(dx_dpsi, psi_vec, initial=0)

        # 绝对深度调整: 底部在 y=-1 (近似，实际应配合 eta 使用)
        y_line = y_from_surf - y_from_surf[-1] - 1.0

        c = np.sqrt(self.F2)
        # 返回: x相对偏移, y绝对深度, u, w (实验室系，向下为负)
        return x_from_surf, y_line, u_arr * c, -v_arr * c



    # =========================================================================
    # 【新增功能】计算任意 (x, sigma) 网格上的速度场
    # =========================================================================
    def calculate_velocity_on_sigma_grid(self, x_req, sigma_req, psi_res=25):
        """
        计算指定 x 和 sigma 网格上的速度矢量。

        参数:
        x_req: 1D array, 目标物理坐标 x 数组
        sigma_req: 1D array, 目标垂向坐标 sigma (0=底部, 1=自由面)
        psi_res: int, 内部计算网格的流函数层数 (越大精度越高，但越慢)

        返回:
        self.data_U: 2D array [len(sigma), len(x)], 水平速度 u
        self.data_W: 2D array [len(sigma), len(x)], 垂直速度 w
        """
        # print(f"正在计算内部全流场 (Phi节点数={len(self.phi_full)}, Psi层数={psi_res})...")

        # 1. 准备计算平面的网格 (Phi, Psi)
        phi_nodes = self.phi_full
        # Psi 从 1 (表面) 向下积分到 0 (底部)
        psi_nodes = np.linspace(1, 0, psi_res)
        d_psi_val = psi_nodes[1] - psi_nodes[0]  # 负值

        # 初始化存储网格 (Complex Plane Grid -> Physical Mapping)
        mesh_shape = (len(psi_nodes), len(phi_nodes))
        X_mesh = np.zeros(mesh_shape)
        Y_mesh = np.zeros(mesh_shape)
        U_mesh = np.zeros(mesh_shape)
        W_mesh = np.zeros(mesh_shape)

        # 2. 计算表面层 (Psi=1) 的物理坐标和速度
        # 重建全域 Tau (利用对称性: tau是偶函数)
        half_len = (len(self.phi_full) + 1) // 2
        tau_half = self.last_tau_right  # 正半轴 tau
        # 拼接得到全域 tau: [tau_N, ..., tau_1, tau_0, tau_1, ..., tau_N]
        # 注意 r_full 的结构: (-r_N... -r_1, 0, r_1... r_N)
        tau_full = np.concatenate((tau_half[1:][::-1], tau_half))

        q_surf = np.exp(tau_full)
        theta_surf = self.theta_full

        # 计算 Jacobian
        jacobian = self.alpha + self.m * np.abs(self.r_full) ** (self.m - 1)

        # 积分表面坐标 (沿 Gamma/Phi 轴)
        dx_surf = (np.cos(theta_surf) / q_surf) * jacobian
        dy_surf = (np.sin(theta_surf) / q_surf) * jacobian

        X_surf = cumulative_trapezoid(dx_surf, self.r_full, initial=0)
        Y_surf = cumulative_trapezoid(dy_surf, self.r_full, initial=0)
        Y_surf = Y_surf - Y_surf[-1]  # 修正基准面

        # 存入第一层
        X_mesh[0, :] = X_surf
        Y_mesh[0, :] = Y_surf
        # 实验室坐标系速度
        U_mesh[0, :] = 1.0 - q_surf * np.cos(theta_surf)
        W_mesh[0, :] = -q_surf * np.sin(theta_surf)  # Tanaka定义theta向下为正，故w = -q*sin

        # 3. 向下积分内部流场 (Layer-by-layer)
        dr = self.r_full[1] - self.r_full[0]

        for i in range(1, len(psi_nodes)):
            psi = psi_nodes[i]
            d_psi = psi - psi_nodes[i - 1]  # 当前层与上一层的差 (负值)

            # --- 3a. 利用解析延拓计算该 Psi 层上的 velocity ---
            # 目标点复势
            W_targets = self.phi_full + 1j * psi
            Phi_source = self.phi_full

            # 向量化计算 Kernel Matrix: Arg_mat[target_idx, source_idx]
            # 形状 broadcasting: (N, 1) - (1, N)
            Arg_mat = np.pi * (W_targets[:, None] - Phi_source[None, :])
            # Kernel = 0.5 * (1 - tanh(x/2))
            Kernel_mat = 0.5 * (1.0 - np.tanh(Arg_mat / 2.0))

            # 积分被积函数: theta(phi') * K * Jac * dr
            Integrand = self.theta_full[None, :] * Kernel_mat * jacobian[None, :]

            # 对源点求和
            Omega = np.sum(Integrand, axis=1) * dr

            tau_layer = Omega.real
            theta_layer = -Omega.imag
            q_layer = np.exp(tau_layer)

            # 存入速度
            U_mesh[i, :] = 1.0 - q_layer * np.cos(theta_layer)
            W_mesh[i, :] = -q_layer * np.sin(theta_layer)

            # --- 3b. 积分物理坐标 (沿 dPsi 方向) ---
            # 公式: dz = (1/q) * exp(i*theta) * i * dPsi
            # dx = - (sin(theta)/q) * dPsi
            # dy =   (cos(theta)/q) * dPsi
            dx_step = -(np.sin(theta_layer) / q_layer) * d_psi
            dy_step = (np.cos(theta_layer) / q_layer) * d_psi

            # 坐标叠加
            X_mesh[i, :] = X_mesh[i - 1, :] + dx_step
            Y_mesh[i, :] = Y_mesh[i - 1, :] + dy_step

        # 4. 插值到用户指定的目标网格 (Interpolation)
        # print("正在插值到目标 (x, sigma) 网格...")

        # 4a. 确定目标点的物理坐标 (Target X, Target Y)
        # 获取 x_req 对应的表面高程 eta
        f_eta = interp1d(X_mesh[0, :], Y_mesh[0, :], kind='linear', fill_value="extrapolate")
        eta_req = f_eta(x_req)

        # 生成 Meshgrid
        XX, SS = np.meshgrid(x_req, sigma_req)
        # 坐标变换: y = -1 + sigma * (1 + eta)
        ZZ = -1.0 + SS * (1.0 + eta_req[None, :])

        # 4b. 执行插值
        # 源数据点云
        points = np.column_stack((X_mesh.flatten(), Y_mesh.flatten()))
        values_U = U_mesh.flatten()
        values_W = W_mesh.flatten()

        # 目标点云
        target_pts = np.column_stack((XX.flatten(), ZZ.flatten()))

        # 使用 griddata (Linear)
        U_interp = griddata(points, values_U, target_pts, method='linear')
        W_interp = griddata(points, values_W, target_pts, method='linear')

        # 重塑形状 [len(sigma), len(x)]
        data_U = U_interp.reshape(XX.shape)
        data_W = W_interp.reshape(XX.shape)
        grid_X = XX
        grid_Z = ZZ

        # print("流场计算完成。")
        return grid_X.transpose(), grid_Z.transpose(), data_U.transpose(), data_W.transpose()

# ================= 主程序入口 =================
def main(solver):
    # 1. 设置工况: 目标相对波高
    eps = solver.eps
    print(f"=== 开始计算: 目标 H/h = {eps} ===")

    # 2. 实例化求解器 (自动执行波高反算)
    # N=80 为高精度计算
    if True:
        # 加载对比数据
        c_exact = np.sqrt(solver.F2)

        # 打印积分物理量
        print("\n--- 积分物理量 (Dimensionless) ---")
        print(f"Froude Number (Fr) : {np.sqrt(solver.F2):.8f}")
        print(f"Mass (M)           : {solver.Mass:.8f}")
        print(f"Impulse (I)        : {solver.Impulse:.8f}")
        print(f"Kinetic Energy (K) : {solver.KineticEn:.8f}")
        print(f"Potential Energy (V): {solver.PotentialEn:.8f}")

        # ================= 绘图 =================
        print("正在绘图...")

    figureLabel = [chr(i) for i in range(ord('a'), ord('z') + 1)]
    if True:

        dirOut = './solitary.Tanaka1986Data'
        if True:

            FigureName = os.path.join(dirOut, f'eps={eps:.2f}')

            goldenRatio = 0.5 * (np.sqrt(5) + 1)

            FigureNx, FigureNy = 3, 2
            FigureNx, FigureNy = 3, 2
            FigureH = 3.2

            fig, axs = plt.subplots(ncols=FigureNy, nrows=FigureNx,
                                    figsize=(FigureH * FigureNy * goldenRatio, FigureH * FigureNx),
                                    sharex=False,
                                    sharey=False,
                                    )

            #
            for FigureIx in range(FigureNx):
                for FigureIy in range(FigureNy):

                    FigureIndex = (FigureIx - 1 + 1) * FigureNy + FigureIy + 1
                    print(f'ix={FigureIx}, iy={FigureIy}, index={FigureIndex}')
                    if FigureNx == FigureNy == 1:
                        myAx = axs
                    elif FigureNx == 1:
                        myAx = axs[FigureIy]
                    elif FigureNy == 1:
                        myAx = axs[FigureIx]
                    else:
                        myAx = axs[FigureIx, FigureIy]

                    myAx.grid()

                    textLabelY = '$YYYY$'
                    textLabelX = '$XXXX$'

                    if FigureIndex == 0:
                        pass
                    elif FigureIndex == 1:
                        # myAx.set_aspect(1.25)
                        textTitle = f'Solitary wave elevation $H/h={eps:.2f}$'
                        textLabelX = '$x/h$'
                        textLabelY = f'$\\eta/h$'

                        tmpData = solver
                        tmp_x = tmpData.grid_x
                        tmp_y = tmpData.data_eta
                        myAx.plot(tmp_x, tmp_y, 'k-', linewidth=2.5, label=f'$z=\\eta$')

                        # tmpData = dc2014
                        # tmp_x = tmpData.x
                        # tmp_y = tmpData.eta
                        # myAx.plot(tmp_x, tmp_y, 'r--', linewidth=2.5, label=f'$z=\\eta$')

                        # Boussinesq 近似解
                        k_bouss = np.sqrt(0.75 * eps)
                        y_bouss = 0.0 + eps * (1.0 / np.cosh(k_bouss * tmp_x)) ** 2
                        myAx.plot(tmp_x, y_bouss, 'b-.', linewidth=1.5, label='Boussinesq')

                        if False:
                            FigureXmin, FigureXmax = np.min(tmp_x), np.max(tmp_x)
                            FigureYmin, FigureYmax = np.min(tmp_y), np.max(tmp_y)
                            FigureTickXN = 10
                            FigureTickYN = 10
                            FigureDx = (FigureXmax - FigureXmin) / FigureTickXN
                            FigureDy = (FigureYmax - FigureYmin) / FigureTickYN
                        else:
                            FigureXmin, FigureXmax = -10, 10
                            FigureXmin, FigureXmax = -2 * solver.L00, 2 * solver.L00
                            FigureYmin, FigureYmax = -1, 1

                            FigureTickXN = 10
                            FigureTickYN = 10
                            FigureDx = (FigureXmax - FigureXmin) / FigureTickXN
                            FigureDy = (FigureYmax - FigureYmin) / FigureTickYN

                        # end
                    elif FigureIndex in [2, 4]:
                        # myAx.set_aspect(1.25)

                        textLabelX = '$x/h$'
                        if FigureIndex == 2:
                            textLabelY = f'$u/c_0$'
                            textTitle = f'Horizontal velocity at $z=-h, \\eta$'

                            tmpData = solver
                            tmp_x = tmpData.grid_x
                            tmp_y = tmpData.data_us
                            myAx.plot(tmp_x, tmp_y, 'k-', linewidth=2.5, label=f'$z=\\eta$' + f', max={np.nanmax(tmp_y):.3f}')

                            tmp_x = tmpData.grid_X[:, -1]
                            tmp_y = tmpData.data_U[:, -1]
                            myAx.plot(tmp_x, tmp_y - 1, 'b-', linewidth=2.5, label=f'$z=-h$' + f', max={np.nanmax(tmp_y):.3f}')

                            # tmpData = dc2014
                            # tmp_x = tmpData.x
                            # tmp_y = tmpData.us
                            # myAx.plot(tmp_x, tmp_y, 'r--', linewidth=2.5, label=f'$z=\\eta$')

                            # tmpData = solver
                            # tmp_x = tmpData.grid_X[:, 0]
                            # tmp_y = tmpData.data_U[:, 0]
                            # print(tmp_x)
                            # myAx.plot(tmp_x, tmp_y, 'm-.', linewidth=2.5, label=f'Tanaka')
                        elif FigureIndex == 4:
                            textLabelY = f'$w/c_0$'
                            textTitle = f'Vertical velocity at $z=\\eta$'

                            tmpData = solver
                            tmp_x = tmpData.grid_x
                            tmp_y = tmpData.data_ws
                            myAx.plot(tmp_x, tmp_y, 'k-', linewidth=2.5, label=f'$z=\\eta$, \nmax={np.nanmax(tmp_y):.3f}')

                            # tmpData = dc2014
                            # tmp_x = tmpData.x
                            # tmp_y = tmpData.ws
                            # myAx.plot(tmp_x, tmp_y, 'r--', linewidth=2.5, label=f'$z=\\eta$')


                        if False:
                            FigureXmin, FigureXmax = np.min(tmp_x), np.max(tmp_x)
                            FigureYmin, FigureYmax = np.min(tmp_y), np.max(tmp_y)
                            FigureTickXN = 10
                            FigureTickYN = 10
                            FigureDx = (FigureXmax - FigureXmin) / FigureTickXN
                            FigureDy = (FigureYmax - FigureYmin) / FigureTickYN
                        else:
                            FigureXmin, FigureXmax = -10, 10
                            FigureXmin, FigureXmax = -2 * solver.L00, 2 * solver.L00
                            FigureYmin, FigureYmax = -1, 1

                            FigureTickXN = 10
                            FigureTickYN = 10
                            FigureDx = (FigureXmax - FigureXmin) / FigureTickXN
                            FigureDy = (FigureYmax - FigureYmin) / FigureTickYN

                        # end
                    elif FigureIndex == 3:
                        # 理论表面极限

                        # myAx.set_aspect(1.25)
                        textTitle = '$c_0=\\sqrt{gh}$'
                        textTitle = f'Horizontal velocity under the wave crest'
                        textLabelX = f'$u/c_0$'
                        textLabelY = f'$z/h$'

                        tmpData = solver
                        tmp_x = tmpData.data_uc
                        tmp_y = tmpData.grid_z
                        myAx.plot(tmp_x, tmp_y, 'k-', linewidth=2.5, label=f'max={np.nanmax(tmp_x):.3f}\nmin={np.nanmin(tmp_x):.3f}')

                        # tmpData = dc2014
                        # tmp_x = tmpData.uc
                        # tmp_y = tmpData.z
                        # myAx.plot(tmp_x, tmp_y, 'r--', linewidth=2.5, label=f'$z=\\eta$')

                        # if hasattr(dataDC, 'uc') and np.any(dataDC.uc):
                        #     myAx.plot(dataDC.uc * np.sqrt(tmpData.F2), dataDC.z, 'b--', linewidth=2, label=f'min={np.nanmin(tmp_x):.3f}, max={np.nanmax(tmp_x):.3f}')

                        if False:
                            FigureXmin, FigureXmax = np.min(tmp_x), np.max(tmp_x)
                            FigureYmin, FigureYmax = np.min(tmp_y), np.max(tmp_y)
                            FigureTickXN = 10
                            FigureTickYN = 10
                            FigureDx = (FigureXmax - FigureXmin) / FigureTickXN
                            FigureDy = (FigureYmax - FigureYmin) / FigureTickYN
                        else:
                            FigureXmin, FigureXmax = 0, 1
                            FigureYmin, FigureYmax = -1, 1

                            FigureTickXN = 10
                            FigureTickYN = 10
                            FigureDx = (FigureXmax - FigureXmin) / FigureTickXN
                            FigureDy = (FigureYmax - FigureYmin) / FigureTickYN

                        # end

                    elif FigureIndex in [ 5, 6]:

                        textLabelX = '$x/h$'
                        textLabelY = f'$z/h$'

                        tmpData = solver
                        tmp_x = tmpData.grid_X
                        tmp_y = tmpData.grid_Z

                        if FigureIndex == 5:
                            tmp_z = tmpData.data_U
                            tmp_cmap = 'viridis'
                            colorBarNumber = 6
                            tmp_level = np.linspace(-0.0, 0.6, 2 * colorBarNumber + 1, endpoint=True)
                            textColorBar = '$u/c_0$'
                            textTitle = f'Horizontal velocity under the wave surface'
                        elif FigureIndex == 6:
                            tmp_z = tmpData.data_W
                            tmp_cmap = 'RdBu_r'
                            colorBarNumber = 4
                            tmp_level = np.linspace(-0.2, 0.2, 2 * colorBarNumber + 1, endpoint=True)
                            textColorBar = '$w/c_0$'
                            textTitle = f'Vertical velocity under the wave surface'

                        # myAx.fill_between(np.arange(-10, 10, 1), -30, 30,  # y下限和上限
                                          # color=gmtColor['land'],
                                          # alpha=1.0,
                                          # label='Impermeable breakwater')

                        CS = myAx.contourf(
                            tmp_x, tmp_y, tmp_z,
                            tmp_level,
                            cmap=tmp_cmap,
                            extend='both'
                        )

                        # myAx.tricontour(tmpData.triang, np.transpose(tmpData.alpha), (0.05,), colors='cyan',
                        #                 linewidths=2)
                        # myAx.tricontour(tmpData.triang, np.transpose(tmpData.alpha), (0.50,), colors='black',
                        #                 linewidths=2)
                        # myAx.tricontour(tmpData.triang, np.transpose(tmpData.alpha), (0.95,), colors='blue',
                        #                 linewidths=2)

                        if False:
                            tmpData = B.exp
                            tmp_x = tmpData[:, 0] - 7
                            tmp_y = tmpData[:, 1] - 7.5

                            myAx.plot(tmp_x, tmp_y, 'ro'
                                      , linewidth=2
                                      , label='Hsiao & Lin (2010)'
                                      , markersize=MarkerSize * 1.5
                                      , markevery=1
                                      )
                            myAx.plot(tmp_x, tmp_y - 999, 'c-',
                                      label='$\\alpha=0.05$',
                                      linewidth=2., markersize=6)
                            myAx.plot(tmp_x, tmp_y - 999, 'k-',
                                      label='$\\alpha=0.50$',
                                      linewidth=2., markersize=6)
                            myAx.plot(tmp_x, tmp_y - 999, 'b-',
                                      label='$\\alpha=0.95$',
                                      linewidth=2.5, markersize=6)

                        # 2. 创建图内嵌入轴 (inset_axes)
                        # [x0, y0, width, height] 均为相对于 myAx 的比例坐标
                        # 如下设置将 colorbar 放在图内右上角
                        # ax_ins = myAx.inset_axes([0.1, 0.9, 0.8, 0.05])
                        # cbar = fig.colorbar(CS,
                        #                     ax=ax_ins,
                        #                     # ticks=[-5, 0, 5, 10, 15, 20, 25],
                        #                     # ticks=tmp_level,
                        #                     ticks=np.linspace(tmp_level.min(), tmp_level.max(), 1 + colorBarNumber,
                        #                                       endpoint=True),
                        #                     # format=mticker.FixedFormatter(['< -1', '0', '> 1']),
                        #                     # shrink=0.9, aspect=20,
                        #                     orientation='horizontal',
                        #                     # location='right',
                        #                     )

                        # 2. 创建图内嵌入轴 (inset_axes)
                        # [x0, y0, width, height] 均为相对于 myAx 的比例坐标
                        # 如下设置将 colorbar 放在图内右上角
                        # ax_ins = myAx.inset_axes([0.5, 0.99, 0.5, 0.03])

                        # 3. 绘制竖直 colorbar
                        cbar = fig.colorbar(CS,
                                            cax=myAx.inset_axes([0.5, 0.95, 0.5, 0.03]),
                                            ticks=np.linspace(tmp_level.min(), tmp_level.max(), 1 + colorBarNumber,
                                                              endpoint=True),
                                            orientation='horizontal',
                                            )

                        cbar.set_label(textColorBar, labelpad=1)

                        # 反转刻度标签顺序
                        cbar.ax.invert_yaxis()  # 这一行仍然需要

                        # if FigureIndex in [1,2]:
                        #     cbar.set_label('$\\rho$/h', labelpad=1)
                        # elif FigureIndex in [3, 4]:
                        #     cbar.set_label('u', labelpad=1)
                        # elif FigureIndex in [5, 6]:
                        #     cbar.set_label('w', labelpad=1)

                        # myAx.text(0.01, 0.99,
                        #           't = %.2f s' % (tmpData.thisTime),
                        #           color='white',
                        #           fontsize='large',
                        #           horizontalalignment='left',
                        #           verticalalignment='top',
                        #           bbox={'facecolor': 'blue', 'alpha': 0.0, 'edgecolor': 'blue', 'pad': 0.5},
                        #           transform=myAx.transAxes,
                        #           )

                        if False:
                            FigureXmin, FigureXmax = np.min(tmp_x), np.max(tmp_x)
                            FigureYmin, FigureYmax = np.min(tmp_y), np.max(tmp_y)
                            FigureTickXN = 10
                            FigureTickYN = 10
                            FigureDx = (FigureXmax - FigureXmin) / FigureTickXN
                            FigureDy = (FigureYmax - FigureYmin) / FigureTickYN
                        else:
                            FigureXmin, FigureXmax = -10, 10
                            FigureXmin, FigureXmax = -2 * solver.L00, 2 * solver.L00
                            FigureYmin, FigureYmax = -1, 1

                        FigureTickXN = 10
                        FigureTickYN = 10
                        FigureDx = (FigureXmax - FigureXmin) / FigureTickXN
                        FigureDy = (FigureYmax - FigureYmin) / FigureTickYN
                    # end

                    # * * * * + * * * * + * * * * + * * * * + * * * * + * * * * + * * * * +
                    if FigureIndex in [1, 2, 3, 4]:
                        myAx.legend(
                            markerscale=1.0,
                            loc="lower right",
                            # fontsize="medium",
                            edgecolor="black",
                            facecolor="white",
                            ncol=1,
                            # title = "Velocity CDF",
                        )

                    if FigureNx * FigureNy > 1:
                        textTitle = f'({figureLabel[FigureIndex - 1]}) ' + textTitle
                        textLabelX += f'\n{textTitle}'
                    # myAx.set_title(textTitle, loc="left")
                    '''
                        WARN: following lines can not be modified, if YOU know how to do them.
                    '''
                    # - - - - + - - - - + - - - - + - - - - + - - - - + - - - - + - - - - +
                    myAx.set_xlim(FigureXmin, FigureXmax)
                    myAx.set_ylim(FigureYmin, FigureYmax)

                    myAx.set_ylabel(textLabelY)
                    myAx.set_xlabel(textLabelX)

                    FigureTickX = np.linspace(FigureXmin, FigureXmax, num=FigureTickXN + 1)
                    FigureTickY = np.linspace(FigureYmin, FigureYmax, num=FigureTickYN + 1)

                    myAx.set_xticks(FigureTickX)
                    myAx.set_yticks(FigureTickY)
                    # - - - - + - - - - + - - - - + - - - - + - - - - + - - - - + - - - - +

            plt.tight_layout()
            plt.savefig(f'{FigureName}.png', bbox_inches='tight', dpi=500)
            plt.savefig(f'{FigureName}.pdf', bbox_inches='tight')
            plt.show()
            # plt.close()


def saveData(solver):


    # 1. 创建保存目录
    output_dir = 'tanaka1986Data'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 2. 保存表面数据 (x, eta, us, ws)
    # 将 1D 数组合并为 2D 矩阵 [N, 4]
    surface_data = np.column_stack((solver.grid_x, solver.data_eta, solver.data_us, solver.data_ws))
    surface_filename = os.path.join(output_dir, f'eps={eps:.2f}surface.csv')
    np.savetxt(surface_filename, surface_data, delimiter=',', fmt='%.6e, %.6e, %.6e, %.6e',
               header='x,eta,us,ws', comments='# surface height, velocity of the solitary wave of Tanaka (1986).\n')
    print(f"已保存表面数据至: {surface_filename}")

    # 3. 保存全流场数据 (grid_X, grid_Z, data_U, data_W)
    # 由于场数据是 2D 矩阵 [N_phi, N_psi]，将其展平为 [TotalPoints, 4] 以便存储
    field_data = np.column_stack((
        solver.grid_X.flatten(),
        solver.grid_Z.flatten(),
        solver.grid_Sigma.flatten(),  # 新增列
        solver.data_U.flatten(),
        solver.data_W.flatten()
    ))
    field_filename = os.path.join(output_dir, f'eps={eps:.2f}Field.csv')
    np.savetxt(field_filename, field_data, delimiter=',', fmt='%.6e, %.6e, %.6e, %.6e, %.6e',
               header='grid_X,grid_Z,grid_Sigma,data_U,data_W',
               comments='# water velocity beneath the surface of the solitary wave of Tanaka (1986).\n')
    print(f"已保存流场数据至: {field_filename}")

    # 3. 追加保存波速汇总数据 (eps, sqrt(F2)) 到 tanakaC.csv
    # 2. 汇总保存波速及所有积分物理量
    eps_val = solver.H
    # 定义汇总文件路径
    prop_path = os.path.join(output_dir, '0EpsCelerityMCIKP.csv')

    # 构造数据行：[振幅, 波速, 质量, 环量, 冲量, 动能, 势能]
    # 对应 MATLAB: SWP(1)到SWP(6) [cite: 312, 1294]
    prop_data = np.array([[
        eps_val,  # eps
        np.sqrt(solver.F2),
        solver.Mass,
        solver.Circulation,
        solver.Impulse,
        solver.KineticEn,
        solver.PotentialEn,
        solver.KineticEn+solver.PotentialEn,
    ]])

    # 检查文件是否存在以决定是否写入表头
    file_exists = os.path.isfile(prop_path)
    with open(prop_path, 'a') as f:
        if not file_exists:
            # 写入表头，方便后续分析
            f.write('eps,c_val,mass,circulation,impulse,kinetic_en,potential_en, All_en\n')

        # 使用 savetxt 追加数据行
        np.savetxt(f, prop_data, delimiter=',',
                   fmt='%.3f, ' + ', '.join(['%.10e'] * 7))

    print(f"--- 物理量汇总已更新 (eps={eps_val:.2f}) ---")
    print(f"路径: {prop_path}")


class classSolitary:
    """
    用于加载 Denys Dutykh (2014) 的高精度数值解数据进行对比。
    """

    def __init__(self, eps):
        base_path = 'dc2014Data'

        try:
            # 1. 加载汇总属性文件 (eps, c_val, mass, circulation, impulse, kinetic_en, potential_en)
            # 对应之前 saveData 中生成的 0EpsCelerityMCIKP.csv
            prop_path = os.path.join(base_path, '0EpsCelerityMCIKP.csv')
            data_prop = np.loadtxt(prop_path, delimiter=',', skiprows=1)  # 跳过表头

            # 找到最接近目标 eps 的行索引
            idx = np.argmin(np.abs(data_prop[:, 0] - eps))

            # 提取无量纲积分量 (假设保存时已经是无量纲或根据需要转换)
            self.c = data_prop[idx, 1]  # 无量纲波速 (Froude数)
            self.Mass = data_prop[idx, 2]
            self.Circulation = data_prop[idx, 3]
            self.Impulse = data_prop[idx, 4]
            self.KineticEn = data_prop[idx, 5]
            self.PotentialEn = data_prop[idx, 6]

            # 2. 加载表面分布数据 (x, eta, u, w)
            # 注意文件名格式需与 saveData 保持一致
            surf_path = os.path.join(base_path, f'eps={eps:.2f}surface.csv')
            data_surf = np.loadtxt(surf_path, delimiter=',', skiprows=1)

            self.x = data_surf[:, 0]
            self.eta = data_surf[:, 1]
            self.us = data_surf[:, 2]
            self.ws = data_surf[:, 3]

            # 3. 加载垂向剖面数据 (假设读取场数据中 x=0 的部分或专门的剖面文件)
            try:
                field_path = os.path.join(base_path, f'eps={eps:.2f}Field.csv')
                data_field = np.loadtxt(field_path, delimiter=',', skiprows=1)

                # 提取 x 接近 0 的点作为垂向剖面 (uc)
                x_tol = 1e-3
                center_mask = np.abs(data_field[:, 0]) < x_tol
                if np.any(center_mask):
                    profile = data_field[center_mask]
                    # 按高度排序
                    sort_idx = np.argsort(profile[:, 1])
                    self.z = profile[sort_idx, 1]
                    self.uc = profile[sort_idx, 3]  # data_U
                else:
                    self.z, self.uc = np.zeros(10), np.zeros(10)
            except:
                self.z, self.uc = np.zeros(10), np.zeros(10)

            print(f"成功加载 DC2014 数据: eps={data_prop[idx, 0]:.3f}, c={self.c:.6f}")

        except Exception as e:
            print(f"Warning: Could not load DC2014 data: {e}")
            # 默认空数据初始化
            self.c, self.Mass, self.Circulation = 1.0, 0.0, 0.0
            self.Impulse, self.KineticEn, self.PotentialEn = 0.0, 0.0, 0.0
            self.x = np.linspace(-10, 10, 100)
            self.eta = np.zeros_like(self.x)
            self.us = np.zeros_like(self.x)
            self.ws = np.zeros_like(self.x)
            self.z = np.linspace(-1, 0, 10)
            self.uc = np.zeros_like(self.z)

if __name__ == "__main__":


    for eps in np.arange(0.79, 0.8, 0.05):
        dc2014 = classSolitary(eps=eps)
        solver = TanakaSolver(eps=eps, N=80)
        # saveData(solver)
        # main(solver)

