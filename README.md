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

确保依赖已安装：

pip install flask rdkit-pypi requests

```bash
git clone https://github.com/Chengkai04/molview.git
cd molview
pip install -r requirements.txt
python app.py
打开浏览器访问：http://127.0.0.1:5001


## 🔄 更新日志

### 2025-07-03

- ✅ 支持通过分子式查询化合物结构（例如 `c10h22`, `c6h6`, `h2o` 等）。
- 🔁 自动规范用户输入的大小写（`c10h22` 会自动变为 `C10H22`）。
- 📡 改进 SMILES 查询逻辑：
  - 优先从名称接口获取
  - 若失败，使用分子式接口获取 CID → 再查 SMILES
  - 支持解析 PubChem 返回的复杂 JSON 数据结构

---

## 🔍 示例输入支持（名称 / 分子式）

| 输入       | 返回结构                         |
|------------|----------------------------------|
| `aspirin`  | ✅ 阿司匹林结构图                 |
| `C10H22`   | ✅ 正癸烷结构图                   |
| `c10h22`   | ✅ 自动规范大小写 → 成功查询       |
| `benzene`  | ✅ 苯的结构图                     |
| `c6h6`     | ✅ 苯（分子式）                   |
| `water`    | ✅ 水                             |
| `h2o`      | ✅ 自动规范大小写 → 成功查询       |

---

## 🚀 部署说明

部署平台建议使用 [Render](https://render.com/) 或本地运行：

```bash
python app.py


2025-07-08
🧬 页面新增展示内容：

CID

分子式（Molecular Formula）

分子量（Molecular Weight）

IUPAC 名称（IUPAC Name）

Canonical SMILES 字符串

📦 新增下载功能：

📥 PNG 图片一键下载

📥 SVG 矢量图一键下载

🎨 页面 UI 大幅美化：

引入 Bootstrap 5 实现响应式、现代化界面

搜索时显示加载提示（loading 提示）

结果区域以卡片方式美观呈现