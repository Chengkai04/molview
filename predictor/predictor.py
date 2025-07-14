import joblib
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import os

# 加载模型
def load_model(model_path="model/solubility_model.pkl"):
    model_path = os.path.abspath(model_path)
    print(f"[DEBUG] 正在加载模型：{model_path}")
    assert os.path.exists(model_path), f"❌ 找不到模型文件：{model_path}"
    model = joblib.load(model_path)
    return model

# SMILES 转换为指纹
def smiles_to_features(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
    return np.array(list(fp))

# 溶解度预测
def predict_solubility(smiles, model):
    features = smiles_to_features(smiles)
    if features is None:
        return None
    prediction = model.predict([features])[0]
    return round(float(prediction), 3)  # 保留三位小数
