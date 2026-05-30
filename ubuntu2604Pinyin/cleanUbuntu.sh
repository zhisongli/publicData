# 1. 清理 apt 缓存（最有效）
sudo apt clean                  # 删除所有已下载的 .deb 包
sudo apt autoclean              # 只删除过期版本的 .deb 包
sudo apt autoremove --purge     # 删除孤立依赖包（含配置文件）

# 2. 清理 journal 日志（保留最近 2 周）
sudo journalctl --vacuum-time=2weeks
sudo journalctl --vacuum-size=200M

# 3. 清理缩略图缓存
rm -rf ~/.cache/thumbnails/*

# 4. 清理 snap 旧版本（如果用 snap）
snap list --all | awk '/disabled/{print $1, $3}' | \
  while read name rev; do sudo snap remove "$name" --revision="$rev"; done


# 5. 查看并清理旧内核（保留当前 + 1 个备用）
# dpkg --list | grep linux-image       # 查看已安装内核
# uname -r                             # 当前运行内核
# sudo apt autoremove --purge          # 通常会自动处理旧内核

# 6. 清理 pip 缓存
pip cache purge
pip3 cache purge

# 7. 清理 conda 缓存（如果用 conda/mamba）
conda clean --all -y

# 8. 清理 npm 缓存
npm cache clean --force 2>/dev/null || true
