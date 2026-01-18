#!/bin/bash

# PDF 转 PNG 脚本
# 将脚本所在目录中的所有 PDF 文件转换为 PNG 格式
# 保持最高分辨率，使用透明背景

set -e  # 遇到错误时退出脚本

# 获取脚本所在目录
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "工作目录: $SCRIPT_DIR"

# 检查是否安装了 ImageMagick
if ! command -v convert &> /dev/null; then
    echo "错误: ImageMagick 未安装！"
    echo "请安装 ImageMagick:"
    echo "  Ubuntu/Debian: sudo apt-get install imagemagick"
    echo "  macOS: brew install imagemagick"
    echo "  CentOS/RHEL: sudo yum install ImageMagick"
    exit 1
fi

# 检查是否安装了 Ghostscript（用于处理某些 PDF）
if ! command -v gs &> /dev/null; then
    echo "警告: Ghostscript 未安装，某些 PDF 可能无法正确处理"
    echo "建议安装:"
    echo "  Ubuntu/Debian: sudo apt-get install ghostscript"
    echo "  macOS: brew install ghostscript"
fi

# 创建输出目录（如果不存在）
OUTPUT_DIR="$SCRIPT_DIR"
if [ ! -d "$OUTPUT_DIR" ]; then
    mkdir -p "$OUTPUT_DIR"
    echo "创建输出目录: $OUTPUT_DIR"
fi

# 查找所有 PDF 文件
pdf_files=("$SCRIPT_DIR"/*.pdf)

echo "找到 ${#pdf_files[@]} 个 PDF 文件:"

# 处理每个 PDF 文件
converted_count=0
failed_files=()

for pdf_file in "${pdf_files[@]}"; do
    # 跳过非文件项（当没有匹配时，通配符会返回自身）
    
    filename=$(basename "$pdf_file" .pdf)
    echo "正在处理: $filename.pdf"
    
    # 为多页 PDF 创建子目录
    if [ 1 = 1 ]; then
        
        output_file="${filename}.png"
        echo "  转换为: $(basename "$output_file")"
        
        # 使用高分辨率设置转换
        convert -density 300 "$pdf_file" \
           -quality 100 \
           -colorspace RGB \
           -background none \
           -alpha remove \
           -alpha off \
           -strip \
           "$output_file"
    fi
    
    echo ""
done
