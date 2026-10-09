# 力學 Project — 牛頓力學的數值積分（Group 22）

古典力學課程的分組專題：以 **velocity Verlet** 與 **RK4** 兩種數值積分方法，模擬**單擺**（小角度近似 vs 完整非線性；無阻尼 vs 一次方阻尼），並比較兩種方法的表現。

- **Part 1**：查詢數值方法、選擇力學問題、設計數值實驗（已完成）
- **Part 2**：依 Part 1 的設計實際執行計算實驗，並撰寫實驗報告（`part2/`）

## 資料夾結構

| 路徑 | 內容 |
|---|---|
| `Project01student.md` / `.pdf` | 作業說明（Part 1 題目） |
| `力學_project1_G22.md` / `.pdf` | 本組 Part 1 報告：方法選擇與四個實驗的設計 |
| `numerical_tools_student_release/student_release/` | 課程提供的數值工具包（**請勿修改**） |
| `part2/` | Part 2 實驗程式、圖表、數值結果與 LaTeX 報告 |
| `part2/README.md` | Part 2 詳細說明：實驗內容、物理模型、與 Part 1 預測的比較、協作規則 |
| `part2/README_for_agents.md` | 給 AI agent 的工作規則（用 AI 協作時請先讓它讀這份） |
| `part2/report/report_G22.tex` | Part 2 實驗報告（繁體中文，XeLaTeX） |

## 快速開始

```bash
git clone git@github.com:glen0414/dynamics_project1.git
cd dynamics_project1
python -m pip install -r numerical_tools_student_release/student_release/requirements.txt  # numpy, matplotlib
cd part2
python run_all.py      # 約 20 秒；圖存到 part2/figures/，數值摘要存到 part2/results/
```

## 編譯報告

```bash
cd part2/report
xelatex report_G22.tex   # 跑兩次；需 xeCJK 與 Noto Serif CJK TC 字型
```

使用 Overleaf 時：上傳 `report_G22.tex`，並把 `part2/figures/` 內的 PNG 放到專案的 `figures/` 資料夾，編譯器選 **XeLaTeX**。

## 主要結論（Part 2）

1. 大角度時非線性週期變長（60° 時 T/T₀ = 1.073），與 Δt 無關，屬於物理現象。
2. 無阻尼時 VV 能量誤差有界、不漂移；RK4 能量單調流失，但在 Δt = 0.01 s 下誤差反而小約 140 倍。
3. 工具包的 velocity Verlet 只接受位置相依的加速度，無法直接處理阻尼；自行修改後，延遲版為一階、隱式版為二階，RK4 為四階。
4. 有阻尼的大角度單擺，頻率會收斂到小角度近似的頻率，但早期累積的相位差不會消失。

詳見 `part2/README.md` 與報告。
