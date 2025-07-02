from flask import Flask, render_template, request, send_file, session
from rdkit import Chem
from rdkit.Chem import Draw
import io
import base64
import requests
import os

app = Flask(__name__)
app.secret_key = "something_secret"

# ✅ 新版查询函数：支持名称和分子式自动判断
def get_smiles(query):
    # 尝试按名称查
    name_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{query}/property/CanonicalSMILES/TXT"
    response = requests.get(name_url)
    if response.status_code == 200:
        smiles = response.text.strip()
        if Chem.MolFromSmiles(smiles):
            return smiles

    # 尝试按分子式查
    formula_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/formula/{query}/property/CanonicalSMILES/TXT"
    response = requests.get(formula_url)
    if response.status_code == 200:
        # 多行结果，每行一个 SMILES
        smiles_list = response.text.strip().split("\n")
        for smiles in smiles_list:
            if Chem.MolFromSmiles(smiles):  # ✅ 能成功解析
                return smiles

    return None



@app.route("/", methods=["GET", "POST"])
def index():
    img_data = None
    error = None

    if request.method == "POST":
        query = request.form["compound"]
        smiles = get_smiles(query)

        if smiles:
            mol = Chem.MolFromSmiles(smiles)
            img = Draw.MolToImage(mol)

            buf = io.BytesIO()
            img.save(buf, format="PNG")
            buf.seek(0)

            img_data = base64.b64encode(buf.getvalue()).decode("utf-8")
            session["image_bytes"] = img_data
        else:
            error = "❌ 没有找到该化合物，请确认名称或分子式是否正确"

    return render_template("index.html", img_data=img_data, error=error)

@app.route("/download")
def download_image():
    if "image_bytes" not in session:
        return "No image available", 400

    image_bytes = base64.b64decode(session["image_bytes"])
    buf = io.BytesIO(image_bytes)
    buf.seek(0)
    return send_file(buf, mimetype="image/png", as_attachment=True, download_name="molecule.png")

if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5000))
    app.run(host="0.0.0.0", port=port)
