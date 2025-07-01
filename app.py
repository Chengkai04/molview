from flask import Flask, render_template, request, send_file, session
from rdkit import Chem
from rdkit.Chem import Draw
import io
import base64
import requests

app = Flask(__name__)
app.secret_key = "something_secret"  # required for session

def get_smiles_from_name(name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/TXT"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text.strip()
    return None

@app.route("/", methods=["GET", "POST"])
def index():
    img_data = None
    error = None

    if request.method == "POST":
        query = request.form["compound"]
        smiles = get_smiles_from_name(query)

        if smiles:
            mol = Chem.MolFromSmiles(smiles)
            img = Draw.MolToImage(mol)

            # 保存图片到内存中
            buf = io.BytesIO()
            img.save(buf, format="PNG")
            buf.seek(0)

            # 编码 base64 用于网页显示
            img_data = base64.b64encode(buf.getvalue()).decode("utf-8")
            session["image_bytes"] = base64.b64encode(buf.getvalue()).decode("utf-8")  # 保存图像内容
        else:
            error = "❌ 没有找到该化合物"

    return render_template("index.html", img_data=img_data, error=error)

@app.route("/download")
def download_image():
    if "image_bytes" not in session:
        return "No image available", 400

    # 从 session 中读取 base64 并转回字节
    image_bytes = base64.b64decode(session["image_bytes"])
    buf = io.BytesIO(image_bytes)
    buf.seek(0)
    return send_file(buf, mimetype="image/png", as_attachment=True, download_name="molecule.png")

if __name__ == "__main__":
    import os
port = int(os.environ.get("PORT", 5000))
app.run(host="0.0.0.0", port=port)

