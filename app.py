from flask import Flask, render_template, request, send_file, session
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
import io
import base64
import requests
import os
import re
from urllib.parse import quote_plus

# 导入溶解度预测相关函数
from predictor.predictor import load_model, predict_solubility


app = Flask(__name__)
app.secret_key = "something_secret"

# 应用启动时加载模型
model = load_model()

COMMON_COMPOUNDS = [
    "water", "ethanol", "glucose", "benzene", "aspirin", "caffeine", "acetone", "sodium chloride"
]

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


def get_compound_info(query):
    """Fetch SMILES first, then retrieve CID / formula / weight.
    Returns dict on success else None."""

    smiles = get_smiles(query)
    if not smiles:
        return None

    # 1️⃣ Get CID from SMILES
    cid = None
    cid_url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/"
        f"{quote_plus(smiles)}/cids/JSON"
    )
    try:
        cid_resp = requests.get(cid_url)
        if cid_resp.status_code == 200:
            cid_list = cid_resp.json()["IdentifierList"].get("CID", [])
            if cid_list:
                cid = cid_list[0]
    except Exception:
        pass

    # 2️⃣ Fetch molecular formula & weight using CID (preferred) else by name
    formula = weight = iupac_name = None
    if cid:
        prop_url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/"
            "MolecularFormula,MolecularWeight,IUPACName/JSON"
        )
        resp = requests.get(prop_url)
        if resp.status_code == 200:
            try:
                props = resp.json()["PropertyTable"]["Properties"][0]
                formula = props.get("MolecularFormula")
                weight = props.get("MolecularWeight")
                iupac_name = props.get("IUPACName")
            except Exception:
                pass

    # Fallback: try name search for formula/weight if still missing
    if (formula is None or weight is None or iupac_name is None):
        query_std = normalize_formula(query)
        prop_url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{query_std}/property/"
            "MolecularFormula,MolecularWeight,IUPACName/JSON"
        )
        resp = requests.get(prop_url)
        if resp.status_code == 200:
            try:
                props = resp.json()["PropertyTable"]["Properties"][0]
                cid = cid or props.get("CID")
                formula = formula or props.get("MolecularFormula")
                weight = weight or props.get("MolecularWeight")
                iupac_name = iupac_name or props.get("IUPACName")
            except Exception:
                pass

    if cid and formula and weight:
        predicted_solubility = predict_solubility(smiles, model)
        return {
            "cid": cid,
            "formula": formula,
            "weight": weight,
            "smiles": smiles,
            "iupac": iupac_name,
            "predicted_solubility": predicted_solubility,
        }

    return None


@app.route("/", methods=["GET", "POST"])
def index():
    img_data = None
    error = None
    info = {}

    # 搜索历史
    history = session.get("search_history", [])

    if request.method == "POST":
        query = request.form["compound"]
        info = get_compound_info(query)
        info = info or {}

        # 追加历史（去重，最新在前，最多10条）
        if query and (not history or query != history[0]):
            history = [query] + [h for h in history if h != query]
            history = history[:10]
            session["search_history"] = history

        if info:
            mol = Chem.MolFromSmiles(info["smiles"])
            img = Draw.MolToImage(mol)
            buf = io.BytesIO()
            img.save(buf, format="PNG")
            buf.seek(0)
            img_data = base64.b64encode(buf.getvalue()).decode("utf-8")
            session["current_smiles"] = info["smiles"]
        else:
            error = "❌ 没有找到该化合物，请确认名称或分子式是否正确"

    return render_template(
        "index.html",
        img_data=img_data,
        error=error,
        cid=info.get("cid"),
        formula=info.get("formula"),
        weight=info.get("weight"),
        canonical_smiles=info.get("smiles"),
        iupac_name=info.get("iupac"),
        predicted_solubility=info.get("predicted_solubility"),
        smiles_q=quote_plus(info.get("smiles") or "") if info else None,
        history=history,
        common_compounds=COMMON_COMPOUNDS,
    )

@app.route("/download")
def download_image():
    smiles = request.args.get("smiles")
    if not smiles:
        smiles = session.get("current_smiles")
    if not smiles:
        return "No image available", 400

    mol = Chem.MolFromSmiles(smiles)
    img = Draw.MolToImage(mol)
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    buf.seek(0)
    return send_file(buf, mimetype="image/png", as_attachment=True, download_name="molecule.png")


@app.route("/download_svg")
def download_svg():
    smiles = request.args.get("smiles")
    if not smiles:
        smiles = session.get("current_smiles")
    if not smiles:
        return "No image available", 400

    mol = Chem.MolFromSmiles(smiles)
    drawer = rdMolDraw2D.MolDraw2DSVG(300, 300)
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    svg_text = drawer.GetDrawingText().encode("utf-8")
    buf = io.BytesIO(svg_text)
    buf.seek(0)
    return send_file(buf, mimetype="image/svg+xml", as_attachment=True, download_name="molecule.svg")

if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5001))
    app.run(host="0.0.0.0", port=port)
