# 🧪 MolView - 分子结构可视化工具

MolView 是一个基于 Flask 构建的 Web 应用，允许用户输入分子名称或化学式，实时查询 PubChem 并生成其结构图像。支持结构图在线查看与下载。

## 🚀 功能特色
- 实时查询分子结构（通过 PubChem API）
- RDKit 可视化分子结构
- 网页展示图像
- 一键下载结构图（PNG）
- 支持中文/英文化合物名或分子式输入

## 🖥️ 技术栈
- Python
- Flask
- RDKit
- HTML + CSS
- requests + base64

## 💡 使用方法

```bash
git clone https://github.com/Chengkai04/molview.git
cd molview
pip install -r requirements.txt
python app.py
打开浏览器访问：http://127.0.0.1:5000