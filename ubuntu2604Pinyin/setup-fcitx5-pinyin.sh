#!/usr/bin/env bash
set -euo pipefail

# ── 1. 安装软件包 ────────────────────────────────────────────────────────────
sudo apt install -y \
  fcitx5 fcitx5-pinyin fcitx5-chinese-addons \
  fcitx5-frontend-gtk3 fcitx5-frontend-gtk4 \
  fcitx5-frontend-qt5 fcitx5-frontend-qt6 \
  fcitx5-config-qt

# ── 2. im-config 选择 fcitx5 ─────────────────────────────────────────────────
im-config -n fcitx5

# ── 3. 环境变量（systemd 用户会话） ──────────────────────────────────────────
mkdir -p ~/.config/environment.d
cat > ~/.config/environment.d/im.conf << 'EOF'
GTK_IM_MODULE=fcitx
QT_IM_MODULE=fcitx
XMODIFIERS=@im=fcitx
EOF

# ── 4. fcitx5 自启动（--replace 替换 ibus） ───────────────────────────────────
mkdir -p ~/.config/autostart
cat > ~/.config/autostart/fcitx5.desktop << 'EOF'
[Desktop Entry]
Name=Fcitx 5
GenericName=Input Method
Comment=Start Fcitx 5 Input Method (replaces ibus)
Exec=fcitx5 --replace -d
Icon=fcitx5
Terminal=false
Type=Application
Categories=System;Utility;
StartupNotify=false
X-GNOME-AutoRestart=false
X-GNOME-Autostart-Notify=false
X-GNOME-Autostart-enabled=true
X-KDE-autostart-after=panel
EOF

# ── 5. 禁用旧版 fcitx4 自启动 ────────────────────────────────────────────────
cat > ~/.config/autostart/fcitx.desktop << 'EOF'
[Desktop Entry]
Hidden=true
EOF

# ── 6. 配置 fcitx5 输入法列表（拼音） ────────────────────────────────────────
mkdir -p ~/.config/fcitx5
cat > ~/.config/fcitx5/profile << 'EOF'
[Groups/0]
Name=Default
Default Layout=us
DefaultIM=pinyin

[Groups/0/Items/0]
Name=keyboard-us
Layout=

[Groups/0/Items/1]
Name=pinyin
Layout=

[GroupOrder]
0=Default
EOF

# ── 7. 立即启动 fcitx5（替换当前会话中的 ibus） ───────────────────────────────
fcitx5 --replace -d
sleep 1

# ── 8. 验证 ──────────────────────────────────────────────────────────────────
echo ""
echo "=== 验证结果 ==="
echo -n "fcitx5 进程: "; pgrep -a fcitx5 || echo "未运行（请重启后再试）"
echo -n "ibus  进程: "; pgrep -a ibus 2>/dev/null || echo "无（正常）"
echo -n "可用输入法: "; fcitx5-remote -l 2>/dev/null || echo "（需重启后生效）"
echo ""
echo "完成。切换输入法快捷键：Ctrl+Space"
echo "如当前应用未生效，注销重新登录后即可。"
