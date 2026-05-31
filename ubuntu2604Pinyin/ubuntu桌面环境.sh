❯ ubuntu2604 安装其他形式的桌面环境后，怎么恢复默认环境的图标？                                                                                                                                                                                                                                                      

● 安装其他桌面环境（如 KDE、XFCE 等）后，Ubuntu GNOME 默认图标可能被替换。以下几种方法可以恢复：

  方法一：通过命令行重置图标主题

  # 恢复 Ubuntu 默认 Yaru 图标主题
  gsettings set org.gnome.desktop.interface icon-theme 'Yaru'

  # 同时恢复 GTK 主题
  gsettings set org.gnome.desktop.interface gtk-theme 'Yaru'
  gsettings set org.gnome.desktop.interface color-scheme 'default'

  方法二：通过 GNOME Tweaks（图形界面）

  sudo apt install gnome-tweaks

  打开 GNOME Tweaks → Appearance → Icons，选择 Yaru。

  方法三：重装 ubuntu-desktop 元包

  sudo apt install --reinstall ubuntu-desktop ubuntu-desktop-minimal

  这会拉回所有默认组件，包括图标包。

  方法四：重置整个 GNOME dconf 配置

  # 重置所有 GNOME 设置（慎用，会重置所有个性化配置）
  dconf reset -f /org/gnome/

  如果桌面图标（回收站、主目录等）不见了

  # 确保 Yaru 图标包已安装
  sudo apt install yaru-theme-icon yaru-theme-gtk

  # 如果用的是 GNOME 扩展管理桌面图标，检查扩展状态
  gnome-extensions list | grep -i desktop
  gnome-extensions enable ding@rastersoft.com  # Desktop Icons NG

  根本原因

  其他 DE（如 KDE）安装时会修改 /etc/alternatives 或 ~/.config 下的默认主题配置，方法一和方法三 通常最快解决问题。

