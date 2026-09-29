# workspace

再生成できる作業用データをまとめます。

- `work/`: Elmer mesh、restart、作業中の中間データ
- `generated/`: project/model から生成する内部 JSON、SIF、図
- `gmsh/`: Gmsh の形状・メッシュ

ここは計算のための作業領域です。保存する解析結果は `../artifacts/`、計算結果と固定実行一式は `../outputs/` に置きます。
