from flask import Flask, render_template, request, send_file, session
from rdkit import Chem
from rdkit.Chem import Draw
import io
import base64
import requests
import os
import re

app = Flask(__name__)
app.secret_key = "something_secret"

def normalize_formula(s):
    return re.sub(r'([a-zA-Z])(\d*)', lambda m: m.group(1).upper() + m.group(2), s)

def extract_smiles_from_pccompound(json_data):
    try:
        props = json_data["PC_Compounds"][0]["props"]
        for prop in props:
            if prop["urn"].get("label") == "SMILES" and prop["urn"].get("name") == "Absolute":
                return prop["value"]["sval"]
    except Exception as e:
        print("❌ 解析 SMILES 失败:", e)
    return None


def get_smiles(query):
    print(f"🧪 [DEBUG] 用户输入：{query}")

    query = normalize_formula(query)
    print(f"🔤 [DEBUG] 标准化后：{query}")

    # 1️⃣ 尝试名称搜索 JSON 结构
    name_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{query}/JSON"
    resp = requests.get(name_url)
    if resp.status_code == 200:
        try:
            json_data = resp.json()
            smiles = extract_smiles_from_pccompound(json_data)
            if smiles and Chem.MolFromSmiles(smiles):
                print("✅ 找到 SMILES:", smiles)
                return smiles
        except Exception as e:
            print("❌ 解析 JSON 失败:", e)

    # 2️⃣ 尝试分子式 → CID → SMILES
    formula_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastformula/{query}/cids/JSON"
    resp = requests.get(formula_url)
    if resp.status_code == 200:
        try:
            cids = resp.json()["IdentifierList"]["CID"]
            if cids:
                cid = cids[0]
                cid_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/JSON"
                cid_resp = requests.get(cid_url)
                if cid_resp.status_code == 200:
                    smiles = extract_smiles_from_pccompound(cid_resp.json())
                    if smiles and Chem.MolFromSmiles(smiles):
                        print("✅ 通过分子式找到 SMILES:", smiles)
                        return smiles
        except Exception as e:
            print("❌ 分子式解析失败:", e)

    print("❌ 没有找到可解析的结构")
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
    port = int(os.environ.get("PORT", 5001))
    app.run(host="0.0.0.0", port=port)
