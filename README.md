# 🧪 MolView - 分子结构可视化工具

MolView 是一个基于 Flask 构建的轻量级 Web 应用，结合 PubChem API 与机器学习模型，支持用户输入分子名称或分子式，实时生成分子结构图并预测其水中溶解度。支持结构图在线展示与一键下载。

🚀 项目特点
输入化合物名称/分子式，自动规范大小写并检索；

后端调用 PubChem API 实时获取分子数据；

用 RDKit 绘制分子结构；

用机器学习模型（RandomForest）预测溶解度；

支持结构图在线查看与下载；

采用 Flask 构建，前端页面美观、简洁、专业。

## 🖥️ 技术栈
-Python
-Flask
-RDKit
-scikit-learn
-PubChem REST API
-HTML + CSS（Bootstrap 5）
-requests + base64

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


📋 更新日志
2025-07-08
🎉 项目重大更新！MolView 2.0 完成：

✅ 集成机器学习模型，实现分子溶解度预测；

✅ 支持输入化合物名称或分子式，获取并显示：

分子结构图

PubChem CID

分子式

分子量

Canonical SMILES

预测溶解度 (log mol/L)

✅ 页面美化为深色科技风，信息以卡片式展示；

✅ 添加副标题，突出机器学习功能；

✅ 优化前后端结构，后端集成已训练好的随机森林回归模型。

