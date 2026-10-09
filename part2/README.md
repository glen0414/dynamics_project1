# 力學 Project Part 2 — 單擺數值實驗（Group 22）

本資料夾是 Part 2（數值計算實驗）的共用工作區。Part 1 報告（`../力學_project1_G22.md`）設計了四個實驗，這裡把它們實際跑出來。

- 數值方法：**velocity Verlet（VV）** 與 **RK4**
- 力學問題：**單擺**（小角度近似 vs 完整非線性；無阻尼 vs 一次方阻尼）
- 數值工具：課程提供的 `../numerical_tools_student_release/student_release/mechanics_integrators.py`（**不修改**）

## 1. 快速開始

```bash
# 需要 numpy、matplotlib（見 toolkit 的 requirements.txt）
cd part2
python run_all.py          # 跑全部實驗，圖存到 figures/，數值摘要存到 results/
python exp2_energy.py      # 也可以單獨跑某一個實驗
```

所有圖都是 PNG（不跳視窗），數值結果是純文字，方便貼進報告。

## 2. 檔案說明

| 檔案 | 內容 |
|---|---|
| `common.py` | 共用參數（g、L、阻尼 γ、Δt）、運動方程式、能量、解析解、週期量測、**自寫的阻尼版 velocity Verlet** |
| `exp1_amplitude.py` | 實驗一：無阻尼，小角度 vs 非線性，θ₀ = 5°、60°（RK4） |
| `exp2_energy.py` | 實驗二：無阻尼非線性、大角度，VV vs RK4 的能量誤差與相空間 |
| `exp3_damped_accuracy.py` | 實驗三：小角度 + 阻尼，VV / RK4 對解析解的誤差；**toolkit 限制的發現** |
| `exp4_damped_nonlinear.py` | 實驗四：有阻尼、大角度，小角度 vs 非線性（RK4） |
| `run_all.py` | 依序執行四個實驗 |
| `figures/` | 輸出的圖 |
| `results/` | 每個實驗的數值摘要（`expN_summary.txt`） |
| `README_for_agents.md` | 給 AI agent 的工作規則（組員用 AI 協作時請讓它先讀這份） |

## 3. 物理模型與符號

角度 θ、角速度 ω = dθ/dt，每單位質量：

- 非線性：θ'' = −(g/L) sin θ − γ ω
- 小角度：θ'' = −(g/L) θ − γ ω
- 能量（每單位質量）：E = ½ L² ω² + gL(1 − cos θ)（小角度模型用 ½ gL θ²）
- 預設參數：g = 9.8 m/s²、L = 1 m、ω₀ = √(g/L)、T₀ = 2π/ω₀ ≈ 2.007 s

## 4. 重要發現：toolkit 的 velocity Verlet 不能直接處理阻尼（實驗三）

toolkit 的 `velocity_verlet(acceleration, t, x0, v0)` 要求 `acceleration(t, x)`，**不能依賴速度**（見 toolkit README 第 10 節）。阻尼 −γω 依賴速度，所以：

1. 直接把 `a(t, x, v)` 丟進 toolkit 的 VV 會出錯（實驗三會實際示範並記錄錯誤訊息）。
2. 因此我們在 `common.py` 自己寫兩個「修改版」VV（不動 toolkit 檔案）：
   - `vv_damped_lagged`：計算 a(n+1) 時用**舊的速度 ωₙ**（最直覺的偷懶寫法）→ 預期只剩一階精度。
   - `vv_damped_implicit`：在速度更新裡把 ω(n+1) 移項解出（阻尼是線性的，可以直接解，不用迭代）→ 保有二階精度。
3. 實驗三比較 RK4、兩種 VV 與解析解，並做 Δt 收斂測試。

報告中可以把這點寫成結論：**「symplectic 的 VV 是為保守力設計的；遇到速度相依的力，必須修改演算法，且修改方式會影響精度。」**

## 5. 如何區分物理與數值誤差（Part 1 第 4 節）

每個實驗都用 Δt 與 Δt/2 各跑一次：

- 結果不隨 Δt 改變 → 物理現象（如大角度週期變長）
- 結果隨 Δt 縮小 → 數值誤差；縮小倍率 ≈ 2^p，p 為方法階數（VV ≈ 4 倍、RK4 ≈ 16 倍）

`results/*.txt` 會直接印出這些比值。

## 6. 目前結果 vs Part 1 預測（Δt = 0.01 s；細節見 `results/`）

| 實驗 | Part 1 預測 | 實際結果 | 是否符合 |
|---|---|---|---|
| 1 | 60° 非線性週期明顯變長 | T/T₀ = 1.073（5° 只有 1.0005），Δt 減半只變 ~1e-8 s，與精確橢圓積分週期一致 | ✅ |
| 2 | RK4 能量流失、相空間內縮；VV 能量守恆 | VV 誤差有界震盪（~2e-4）；RK4 單調流失但**小得多**（2000 s 後 1.3e-6）。內縮只有在粗 Δt = 0.1 才看得到 | ⚠️ 定性對，定量上 RK4 反而更好 |
| 2 | Δt 減半 → RK4 能量流失 ×1/16 | 實測 ×1/32（流失率 ∝ h⁵）；VV ×1/4 | ⚠️ 修正為 1/32 |
| 3 | VV 因速度延遲有累積殘差；RK4 殘差極小 | toolkit VV 直接報錯；延遲版 VV 一階（p = 1.07）、隱式版 VV 二階（p = 2.00）、RK4 四階（p = 4.00） | ✅ 並補充 toolkit 限制 |
| 4 | 非線性頻率收斂到線性頻率，兩者同步停下 | 週期由 2.22 s 收斂到 T_d = 2.009 s ✅；但前期累積的相位差（~0.8 s）**不會消失** | ⚠️ 頻率同步、相位不同步 |

## 7. 協作規則

- **不要修改** toolkit 的 `mechanics_integrators.py`；需要新方法就寫在 `common.py`。
- 參數只在 `common.py` 改，讓所有實驗一致；改了請在 commit 訊息註明。
- 圖的文字用英文（避免中文字型缺失）；報告文字再用中文說明。
- 改完程式請重跑 `python run_all.py`，確認圖與摘要都有更新。
