# Streamlit Conformer/QM Automation App Requirements

## 0. ドキュメントの目的

本ドキュメントは、以下の Python / Streamlit アプリケーションの要件定義書である。

目的

有機分子（SMILES）を入力として、

- コンフォマー生成（RDKit）
- MM 最適化
- UMA ポテンシャルによる構造最適化（オプション）
- GAMESS による QM 計算
- QM + 幾何情報の統合ポスト処理
- DB 更新用バッチの生成

を自動で行う。

スコープ

- GUI: Streamlit アプリ（app.py）
- ライブラリ: RDKit, fairchem, ASE, GAMESS（外部ソフト）, pandas, numpy, etc.
- 付帯ツール: QM リカバリなどの CLI スクリプト（tool/ 配下）

本書を見れば、別の開発者 / 生成AI が同等機能を持つアプリケーションを実装できることを目標とする。

## 1. 想定環境

- OS: Windows（GAMESS が C:\Users\Public\gamess-64 等にインストールされている前提）
- Python: 3.9–3.11 程度

外部プログラム

- GAMESS 実行用の run.bat と gms_progress.ps1
- GAMESS scratch はグローバル（例: C:\Users\Public\gamess-64\scratch）

Python ライブラリ

- streamlit
- pandas
- numpy
- rdkit
- fairchem + FAIRChemCalculator
- ase（ase.io, ase.optimize.LBFGS）

実行形態

- streamlit run app.py で起動

## 2. ディレクトリ構成（論理）

プロジェクトルート（以下 ROOT と書く）は以下のような構造を持つ。

```
ROOT/
├─ app.py                 # Streamlit アプリ本体
├─ output/                # 計算途中の全中間ファイル（ID ごとにフォルダ）
├─ result/                # 解析・要約結果（ID ごとフォルダ）
├─ template/
│   ├─ gms.inp            # GAMESS テンプレート（%COMMENT_SMILES_N% 等のプレースホルダ付き）
│   ├─ run.bat            # GAMESS 起動用バッチ
│   └─ gms_progress.ps1   # 進捗表示用 PowerShell スクリプト
├─ conformer_app/
│   ├─ __init__.py
│   ├─ ui.py              # Streamlit UI 構築
│   ├─ core/
│   │   ├─ config.py      # AppConfig / EmbedConfig / MMConfig / UMAConfig / QMConfig / OutputConfig など
│   │   ├─ runner.py
│   │   ├─ status.py      # StageState, JobProgressRow 等
│   ├─ app_io/
│   │   ├─ excel_io.py    # progress.xlsx の入出力
│   │   ├─ paths.py
│   │   ├─ resume.py
│   ├─ mm_layer/          # コンフォマー生成 + MM
│   ├─ uma_layer/         # UMA 最適化
│   ├─ qm_layer/          # GAMESS 入力生成・実行・解析
│   ├─ extract/           # Geom + QM 統合解析（Geom_QM_summary）
│   └─ rdkit/             # フラグメント集計・RDKit記述子・DB更新ロジック
└─ tool/
    ├─ recover_geom_qm.py # Geom_QM_summary の再生成用ツール
    └─ recover_qm_calc.py # ある ID のみ QM 再実行するリカバリツール
```

実装時はこの構成を基本とするが、厳密に同一である必要はない。ただし、
output/<ID> / result/<ID> / progress.xlsx の概念は必須。

## 3. データフロー概要

- ユーザーが Excel ファイル（requirement.xlsx または progress.xlsx）をアップロード。
- ID と SMILES の一覧を抽出。
- ユーザーが GUI で各種設定（Embed/MM/UMA/QM/Output）を行う。
- 「計算実行」ボタン押下で、選択された ID 群に対し以下を実行（ID ごと直列）。

1. コンフォマー生成（Embed）
2. MM 最適化 & エネルギー / RMSD による構造選抜 → post_MM/
3. UMA ON の場合のみ：UMA 最適化 → UMA/ & UMA_summary.xlsx
4. QM（GAMESS）入力生成 → 実行 → 解析 → QM_summary.xlsx / result/<ID>/<ID>_best.xyz
5. Geom + QM 統合（Geom_QM_summary.xlsx） → result/<ID>/
6. DB 更新用バッチの更新

各 ID ごとの進捗は progress.xlsx に逐次反映される。

run_time.log / total_run_time.log に時間情報を追記。

## 4. 入力 Excel の仕様

### 4.1 requirement.xlsx（新規計算用）

想定フォーマット

- A列: ID
- B列: SMILES
- C列以降: 無視

ヘッダー行の有無は任意だが、B列に SMILES という文字列を含む行はスキップ。

ID が空の場合は、その行の SMILES を ID 代わりにする。

アプリ側の動作:

- parse_requirement_excel(uploaded) にて ID / SMILES を抽出。
- progress.xlsx を新規に作成し、各 ID に対するステージステータス（MM/UMA/QM）や Formal_Charge を保持できるようにする。
- requirement 読み込み後に、以下を自動実行:
  - RDKit フラグメント集計 → result/ 配下に出力
  - RDKit 記述子集計 → result/ 配下に出力（これらは DB 更新用の下準備）

### 4.2 progress.xlsx（途中再開用）

既に一度計算した結果を含む Excel。

必須カラム:

- ID
- SMILES
- MM_status
- UMA_status
- QM_status

Formal_Charge 列があれば、以降の UMA/QM 電荷設定のソースとして利用。

アプリ側の動作:

- アップロードされたファイルに上記 status カラムが存在する場合、
  - それを ./progress.xlsx として保存（上書き）。
- 1列目・2列目から ID / SMILES を抽出。
- progress.xlsx を読み込み、
  - MM_status, UMA_status, QM_status が すべて "DONE"（大文字・小文字無視）となっている ID を 「完了済み」 と判定し、今回の実行対象から除外する。
  - 1つでも DONE でないステージがある ID は、Embed から再実行してよい（現仕様では Embed ～ QM まで通しで実行）。

## 5. Streamlit UI 要件

### 5.1 メイン画面構成

タイトル

- 「DB更新用 自動計算アプリ」
- 「～MM / UMA / QM 自動計算 + DBバッチ更新～」

サイドバー

- 進捗バー（0.0–1.0）
- 即時中断ボタン
  - 押されたら ABORT_FLAG ファイルを作成。
  - UMA / QM 実行中はこのフラグを頻繁にポーリングして即時中断する。
  - フラグ状態を表示（ON / OFF）

本文

- セクション 1: 「Excel 入力」
  - ファイルアップローダ（requirement.xlsx または progress.xlsx）
  - 読み込んだ件数やモード（新規 or resume）を表示
- セクション 2: 「設定」
  - Embed 設定
  - MM 設定
  - UMA 設定（＋UMA使用 ON/OFF）
  - QM 設定
  - 出力設定
- セクション 3: 「実行」
  - 実行ボタン
  - 実行方針の簡単な説明（UMA ON: MM→UMA→QM、OFF: MM→QM など）
  - 入力プレビュー（ID / SMILES）
  - 結果テーブル（ID ごとのステータス / エラーメッセージ）

### 5.2 Embed 設定 UI（コンフォマー生成）

EmbedConfig のフィールド例:

- num_confs: 最大コンフォマー数（初期値 1000）
- prune_rms: 生成時重複除去 RMSD (Å, 負値で OFF)
- etkdg_version: "ETKDGv3" / "ETKDGv2" / "ETKDG"
- random_seed: 乱数シード（-1 でランダム）
- small_ring: Small ring torsions を使うか
- macrocycle: Macrocycle torsions を使うか
- frag_gap_min, frag_gap_max: '.' で分割される複数フラグメント間の表面距離レンジ [Å]

UI は該当する項目を number_input / selectbox / checkbox で提供。

### 5.3 MM 設定 UI

MMConfig の想定フィールド:

- ff_mode: 力場 ("MMFF94s->UFF" / "MMFF94" / "UFF")
- max_iters: MM 最適化反復数

エネルギー重複除去設定:

- equal_tol: エネルギー同値判定 ΔE (kcal/mol)
- energy_window: Emin + window まで許容するエネルギー窓 (kcal/mol)
- rmsd_min: （実装側で使用）構造重複判定 RMSD 閾値 (Å)
- max_keep: UMA に渡す最大構造数（UMA ON 時の上限。デフォルト 10）

要求:

- UMA ON 時:
  - MM 後の構造重複除去（energy + RMSD）後、エネルギー昇順で max_keep 個を UMA に渡す。
- UMA OFF 時:
  - 同様に重複除去したうえで、post_MM に 最大 3 構造のみ を出力する。

### 5.4 UMA 設定 UI

UMAConfig フィールド:

- model_name or model_tag: fairchem のモデル名（例: "uma-s-1p1"）
- device: "cpu" / "cuda"
- fmax: LBFGS 収束閾値 (eV/Å)
- max_steps: LBFGS 最大ステップ数
- charge: 系の総電荷
- spin: スピン多重度 (2S+1)
- task_name: "omol" など
- debug: True/False

UI要件:

- UMA 設定 expander 内に 「UMA を使用する」チェックボックス を設ける。
- ui_uma_block() は (UMAConfig, use_uma: bool) のタプルを返す。

アプリ要件:

- Formal_Charge が progress.xlsx に存在する場合は、
  - UMAConfig.charge を その値で上書きして UMA を実行する。

### 5.5 QM 設定 UI

QMConfig フィールド:

- template_path: gms.inp テンプレートパス（空なら template/gms.inp を探索）
- run_bat_path: GAMESS 実行用 bat パス（空なら template/run.bat などを探索）
- ps1_path: gms_progress.ps1 パス
- ncpus: 1 ジョブあたりの CPU 数
- max_parallel: 将来用。現仕様では使用しないがフィールドは持つ。
- icharge: 電荷（フォールバック）。Formal_Charge が取れない場合に使う。
- mult: 多重度
- pointgroup: ポイントグループ
- smiles_comment: テンプレート中 %COMMENT_SMILES_N% に展開する文字列

アプリ要件:

- Formal_Charge が存在する場合は、
  - QM 実行時には QMConfig.icharge を その値で上書き して用いる。

### 5.6 出力設定 UI

OutputConfig フィールド:

- base_dir: 出力ベースフォルダ（デフォルト ./output）
- move_to_com_origin: XYZ 出力時に重心を原点へ移動するか

## 6. 計算フロー（ID ごと）

以下は 1 つの ID に対する処理フローである。

### 6.1 事前準備 & 即時中断チェック

abort_now() は ABORT_FLAG (./abort_now.flag) の存在で判定。

ID 処理の最初に abort_now() をチェックし、ON の場合は

- MM/UMA/QM のステータスを ERROR にしてスキップ。
- run_time.log にその旨 note を記録。

各 ID ごとに

- t_start: 処理開始時刻を記録

後述する各ステージごとに t_mm_start などを記録し、最後に run_time.log を出力。

### 6.2 Formal Charge の決定

progress.xlsx に Formal_Charge 列が存在する場合:

- ID 行から該当値を取得し、int として解釈。
- 取得成功時は、UMA および QM における電荷設定を これで上書き。

取得できない場合:

- UMAConfig.charge / QMConfig.icharge に設定された値を使用。

### 6.3 MM ステージ

実装関数: run_mm_for_job(mol_id, smiles, embed_cfg, mm_cfg, output_cfg)（既存仕様を踏襲）

処理:

- RDKit によりコンフォマー生成（EmbedConfig）
- MM 最適化（MMFF/UFF 等）
- エネルギー + RMSD による重複除去
- エネルギー窓 (energy_window) によるフィルタリング
- UMA 用 / QM 用候補構造の出力（output/<ID>/MM や post_MM）

戻り値 MMResult（想定フィールド）

- n_generated: 生成コンフォマー数
- n_optimized: MM 最適化に成功した数
- n_candidates: 重複除去後の候補数
- n_selected: 最終的に post_MM に出力した数
- xyz_files: UMA に渡すべき .xyz ファイルリスト（post_MM のパス）

UMA ON / OFF による挙動差

- UMA ON:
  - post_MM には、mm_cfg.max_keep（デフォ 10）までの候補を出力。
- UMA OFF:
  - post_MM には 最大 3 構造 のみ出力（エネルギー昇順＆重複除去済み）。

即時中断:

- MM ステージ内では粒度の細かい abort チェックは必須ではない（現仕様）。
- ID 間やステージ間で abort_now() をチェックし、中断後の ID には入らないようにする。

### 6.4 UMA ステージ（任意）

#### 6.4.1 UMA 実行条件

use_uma = True の場合のみ UMA を実行。

それ以外の場合：

- UMA をスキップし、row.uma_status = DONE とする。
- 後段 QM では UMA を使わず post_MM から構造を取る。

#### 6.4.2 UMA 入力構造の決定

UMA 実行関数：

```
run_uma_for_job(
    mol_id: str,
    mm_result: Any,
    uma_cfg: UMAConfig,
    out_cfg: Optional[OutputConfig],
    is_abort_now: Optional[Callable[[], bool]],
) -> UMAResult
```

UMA に渡す xyz は次の優先順位で決定:

- mm_result.xyz_files があればそれを使用（post_MM/*.xyz を想定）
- 無ければ output/<ID>/post_MM/*.xyz
- それでも無ければ output/<ID>/MM/*.xyz（バックアップ用）

#### 6.4.3 UMA 内部挙動

UMA は fairchem + ASE を用いて以下を実行:

各 xyz について：

- ASE で atoms を読み込み
- atoms.info["charge"] / atoms.info["spin"] を UMAConfig に従って設定
- FAIRChemCalculator を atoms.calc に設定
- LBFGS による最適化（1ステップごとに energy & forces を取得・ログ出力）
- f_max < fmax_thresh で収束判定。max_steps に達したら打ち切り。

出力:

- output/<ID>/UMA/ に UMA 最適化後 xyz を保存（元ファイル名を引き継ぐ）
- output/<ID>/logs/UMA_*.log にステップごとの E/f_max/f_rms などを保存
- UMA_summary.xlsx（ID 直下）を作成

カラム例:

- idx
- xyz_in
- xyz_out
- status: "ok" / "error" / "aborted"
- energy_init_eV
- energy_final_eV
- rmsd_final_A
- message

UMAResult:

- status: "ok" / "error" / "aborted"
- message: エラーメッセージなど
- uma_summary_path: UMA_summary.xlsx のフルパス
- uma_dir: UMA 出力ディレクトリ
- n_structures: UMA 対象構造数

即時中断:

- is_abort_now を LBFGS ループ内で毎ステップ確認。
- フラグ ON の場合は KeyboardInterrupt 相当の処理を行い、
  - 当該構造のステータスを "aborted" として結果に残す。

### 6.5 QM ステージ（GAMESS）

#### 6.5.1 QM 入力構造の選抜

QM 実行関数：

```
run_qm_for_job(
    mol_id: str,
    uma_summary: Union[str, UMAResult, None],
    qm_cfg: Any,
    out_cfg: Any,
    is_abort_now: Optional[Callable[[], bool]],
) -> QMResult
```

UMA ありの場合:

- uma_summary は UMAResult か UMA_summary.xlsx パス。
- UMA_summary.xlsx を読んで、status == "ok" かつ xyz_out が存在する行を抽出。
- energy_final_eV で昇順ソート。
- RMSD クラスタリング:
  - 最安定構造を代表構造として選択。
  - 残りの候補について、既に選ばれた代表構造群との RMSD を計算。
  - RMSD < 0.4 Å のものは同一クラスと見なしてスキップ。
  - 0.4 Å 以上であれば新しいクラスとして代表に採用。
- クラスタ代表構造の中からエネルギー昇順で 最大 3 構造 を選び、
  - それぞれを QM/gms1, QM/gms2, ... 配下に geom.xyz としてコピー。

UMA なしの場合:

- uma_summary is None として呼び出し。
- output/<ID>/post_MM/*.xyz から 最大 3 構造 を N 個の job として QM/gmsN/geom.xyz にコピー。
  - （post_MM 側で既にエネルギー+RMSD 重複除去済みである前提）

#### 6.5.2 GAMESS 入力ファイル生成

テンプレート gms.inp には %COMMENT_SMILES_N% 等のプレースホルダを持たせる。

make_inp_for_slots(id_dir, gms_template_path, smiles_comment, pointgroup, icharge, mult, slots) を用いて、
QM/gms1/gms1.inp のようにスロット名に対応する .inp を生成する。

#### 6.5.3 GAMESS 実行 + NSTEP リスタート + 即時中断

実行は _run_gamess_jobs_sequential(id_dir, qm_cfg, slots, is_abort_now) に委譲。

1 スロットごとに以下を行う：

即時中断チェック

- スロット開始前に is_abort_now() を確認。
- ON の場合は RuntimeError("GAMESS jobs aborted...") を投げて中断。

プロセス起動

- cmd /c run.bat <inp> <ncpus> <logs_dir> <ps1> で実行。
- subprocess.Popen で起動し、communicate(timeout=1.0) をループ。
- 1 秒ごとに is_abort_now() を確認。

ON の場合:

- proc.terminate() を試みる。
- Windows では taskkill /PID <pid> /T /F で子プロセスごと kill を試みる。
- qm_runner.log に「aborted by user」を記録し、RuntimeError を投げる。

returncode / .out チェック

- returncode != 0 → RuntimeError。
- .out ファイルが存在しない → RuntimeError。

収束判定

- _parse_equilibrium_geometry_and_energy(out_file) で
  - TOTAL ENERGY
  - EQUILIBRIUM GEOMETRY LOCATED ブロックの座標
  を取得。
- 座標が得られれば収束とみなし、次の slot へ。

NSTEP 打ち切り判定 & リスタート

- EQUILIBRIUM GEOMETRY が見つからない場合、
  - _check_gamess_max_step(out_file) でメッセージ判定。
  - NSTEP（最大ステップ到達）に起因する終了で ない 場合 → 即 RuntimeError。
- NSTEP が原因のときのみリスタート対象。

リスタート:

- _parse_last_geometry(out_file) で最後の COORDINATES OF ALL ATOMS ARE (ANGS) ブロックから座標を抽出。
- それを QM/gmsN/geom.xyz として上書き。
- make_inp_for_slots(..., slots=[slot]) で当該 slot の gms.inp を再生成。
- 最大 max_restarts 回（例: 2 回）リスタートを試みる。
- それでも収束しなければ RuntimeError。

スクラッチファイルの扱い（要件）

各 GAMESS 実行の終了後（正常/異常に関わらず）、
「同じジョブ名に対応する scratch ファイルが後続ジョブの計算を乱さないようにする」こと。

実装方法の一例:

- run.bat 側で JOB 名ごとのサブディレクトリを使い、終了後にそのディレクトリを削除する。
- もしくは Python 側で slot ごとにユニークな scratch パスを設定し、終了後にそのディレクトリを削除。

目的:
forrtl: severe (24): end-of-file during read, unit 22, file ... gms1.F22.001
のような 異なる計算の scratch 混在エラーを防止 する。

#### 6.5.4 QM 結果のまとめ

_postprocess_qm_for_id(mol_id, id_dir, slots, result_root) で以下を行う：

各 slot の .out から

- 最終エネルギー（TOTAL ENERGY）
- EQUILIBRIUM GEOMETRY の座標

を抽出し、gmsN_final.xyz を出力。

QM_summary.xlsx を ID 直下に生成

カラム例: slot, out_path, energy_hartree, final_xyz

最もエネルギーの低い構造の final_xyz を

result/<ID>/<ID>_best.xyz としてコピー。

戻り値 QMResult フィールド:

- id
- n_jobs
- summary_path (QM_summary.xlsx)
- best_energy_h
- best_xyz_path (result/<ID>/<ID>_best.xyz)

## 7. Geom + QM 統合ポスト処理

目的: QM 計算結果と分子形状・表面情報を統合した
Geom_QM_summary.xlsx を生成する。

### 7.1 対象フォルダの決定

find_qm_target_dir(mol_id) により、

- output/<ID>/QM/gms* のうち
  - .out と *_final.xyz の両方が存在するサブフォルダを候補とし、
  - 更新日時が最も新しいものを 1 つ選ぶ。

見つからなければ、

- 「Geom+QM 統合をスキップしました」という警告を出し、
- note を run_time.log にも記録。

### 7.2 Geom + QM 統合ロジック

関数: run_geom_qm_extract(qm_target_dir: Path) -> Optional[Path]

要件（一部は簡略化した仕様）:

qm_target_dir 配下の *_final.xyz と .out をペアとして処理。

geom_props_logic.analyze_xyz_file(xyz_path) を用いて、

- 原子記号・座標
- SASA（総表面積 / 極性表面積）
- vdW 表面積・体積
- 慣性半径 Rg
- 卵形度 Ovality

などを計算。

元の .out から QM 情報（エネルギー等）を抽出。

これらを 1 行にまとめ、全構造分を表にして Geom_QM_summary.xlsx として出力。

### 7.3 XYZ パーサの要件（大文字小文字）

analyze_xyz_file 内部の read_xyz は、

- 原子記号の大文字・小文字を問わず受け入れる。

例: "si", "Si", "SI" → 元素「Si」として扱う。

RDKit の周期表に渡す前に、

- 1 文字元素: 1文字目を大文字
- 2 文字元素: 1文字目大文字＋2文字目小文字

に正規化すること。

これにより、Element 'SI' not found のようなエラーを回避する。

### 7.4 出力のコピー

run_geom_qm_extract が返した Geom_QM_summary.xlsx を、

result/<ID>/Geom_QM_summary.xlsx としてコピーする。

## 8. DB 更新用バッチ

rdkit.db_update_logic.update_batches_with_geom(result_root, export_root) を呼び出し、以下を行う：

- result/ 配下の Geom_QM_summary 等を元に、
  - db_update_work/ のような作業領域にバッチファイル（Excel 等）を作成・更新。
- 一定数のバッチが溜まったら、
  - ネットワークドライブ DB_EXPORT_ROOT に自動コピーする。

※ 詳細なフォーマットは既存ロジックに依存するため、ここでは概念レベルに留める。

## 9. resume 機能の仕様

起動時／Excel 読み込み時に progress.xlsx をロード。

- MM_status, UMA_status, QM_status がすべて "DONE" の ID は今回の run から除外。
- それ以外の ID は、Embed/MM/UMA/QM を 最初から再度実行 してよい（既存の output/ は上書き可）。

各 ID 処理後、write_progress_excel(rows, "./progress.xlsx") により、

- 対象 ID の行を書き換え（ステータス / エラー / Formal_Charge など）。

実行中にアプリが落ちたり、UMA/QM が途中でこけた場合でも、

- progress.xlsx をアップロードすれば、その状態から再開できる。

## 10. ログ・時間計測

### 10.1 ID ごとの run_time.log

場所: result/<ID>/run_time.log

内容:

```
=== ID <ID> ===
Start_time = YYYY-MM-DD HH:MM:SS
total_time = HH:MM:SS (xxxx.xxx sec)
MM_time = ...
UMA_time = ...
QM_time = ...
note = ... （任意。エラー内容やリカバリ情報など）
```

時間フォーマット:

- _fmt_hms_and_sec(dt) を用いて、
  - datetime.timedelta 由来の H:MM:SS 表記と
  - (1234.567 sec) の両方を併記。

total_time は、

- MM/UMA/QM のうち 最後に終了したステージの終了時刻 - t_start で計算。

note の例:

- "UMA disabled for this job."
- "Geom+QM 統合ポスト処理でエラー: ..."
- "QM recovered from previous failure." など

### 10.2 全体の total_run_time.log

場所: ROOT/total_run_time.log

各実行（run ボタン押下）につき 1 行追記:

```
Start: <start_str> End: <end_str> Elapsed: HH:MM:SS (xxxx.xxx sec) Jobs: <n_jobs>
```

## 11. 即時中断（Abort）仕様

フラグファイル: ./abort_now.flag

サイドバーの「即時中断」ボタンで作成。

UMA / QM レイヤーに is_abort_now コールバックを渡し、

- UMA: LBFGS ステップごとにチェック
- QM: GAMESS 実行中に 1 秒ごとにチェック（communicate(timeout=1)）

中断時の挙動:

UMA:

- 現在処理中の構造を "aborted" と記録し、それ以上の構造は計算しない。
- UMAResult.status は "aborted" or "ok" + 途中まで のようになる。

QM:

- 現在実行中の GAMESS プロセスを terminate + taskkill で kill。
- その slot は RuntimeError として上層へ伝播。
- 以降の slot / ID には進まない。

app.py 側:

- 例外を catch したら、その時点でのステージを ERROR としてマーク。
- run_time.log には「即時中断フラグにより停止」と note を書く。
- progress.xlsx にも ERROR として記録。

## 12. リカバリ用ツール（CLI）

### 12.1 recover_geom_qm.py

目的:

何らかの理由で Geom_QM_summary.xlsx が生成されなかった ID について、

既に存在する QM 出力 (output/<ID>/QM/... や result/<ID>/<ID>_best.xyz) を利用して

Geom + QM 統合を 後から再実行 する。

仕様（例）:

- 引数: python tool/recover_geom_qm.py <path_to_best_xyz>
- 例: python tool/recover_geom_qm.py result/4/4_best.xyz

best.xyz のコメント行（2行目）から

- ID, slot 名（例: 4 gms1 E=...）を読み取る。

output/<ID>/QM/<slot>/ を特定し、その中の .out & *_final.xyz を利用して、

run_geom_qm_extract() 相当の処理を実行。

結果の Geom_QM_summary.xlsx を result/<ID>/ に出力（上書き可）。

### 12.2 recover_qm_calc.py

目的:

途中で QM 計算がこけた ID について、

既に存在する output/<ID>/QM/gms*/gms*.inp を使って QM-only を再実行し、

QM_summary / Geom_QM_summary / result/<ID>/ を再構成する。

仕様（例）:

- 引数: python tool/recover_qm_calc.py <ID>
- 例: python tool/recover_qm_calc.py 8

処理:

- output/<ID>/QM/gms*/ を走査し、gms*.inp が存在する job を対象とする。
- 通常の QM 実行と同様に GAMESS を実行（NSTEP リスタート / 即時中断チェックも同じ）。
- _postprocess_qm_for_id() を呼び出して QM_summary / <ID>_best.xyz を生成。
- run_geom_qm_extract() を呼び Geom_QM_summary を再生成。
- result/<ID>/run_time.log に「QM recovered」的な note を追加。

実装上の注意:

- conformer_app を import できるように、sys.path に ROOT を追加する。
- UMA など他ステージへの依存は極力避け、QM とポスト処理に限定する。

## 13. 非機能要件

再現性

- 必要に応じて random_seed を固定すれば同じコンフォマー集合を再生成可能。

拡張性

- 今後 QM の並列化（QMConfig.max_parallel を利用）や別 QM ソフトへの拡張が可能なよう、
  - QM レイヤーは run_qm_for_job / _run_gamess_jobs_sequential / _postprocess_qm_for_id の 3 層構造に分ける。

ロバスト性

- 各ステージで例外発生時には、
  - ステータスを ERROR に設定
  - エラーメッセージを error_message に保存
  - 可能な限り他の ID の処理は継続

ログ

- logs/ ディレクトリ配下のログ（qm_runner.log, UMA_xxx.log など）を詳細に残し、
  - 後から問題解析できるようにする。
