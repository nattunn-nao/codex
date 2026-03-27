# CMP吸着構造探索アプリ

Streamlit GUI で、事前登録済み Si/SiO2 スラブに対して添加剤分子の吸着構造を探索するアプリです。

## 主な機能
- 表面モデルのレジストリ選択
- `adsorbates/` 内 `.xyz` の逐次処理
- 配置グリッド(平行移動・高さ・回転)生成
- 幾何フィルタ
- GFN-FF 粗最適化 + 枝刈り
- UMA(OC25) 本最適化
- ランキング表示/CSV・JSON保存

## セットアップ
```bash
cd cmp_adsorption_app
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

## 起動
```bash
bash scripts/run_app.sh
```

## CLI実行
```bash
python -m cmp_adsorption.cli --surface si001
```

## 注記
- GFN-FF(xTB)・UMA(OC25) 実体がない環境では、`mock=true` で擬似エネルギー評価にフォールバックします。
- `data/surfaces/<name>/fixed_index.txt` または `fixed_indices.txt` を自動認識します。
