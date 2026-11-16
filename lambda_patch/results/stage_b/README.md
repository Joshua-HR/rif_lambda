# 단계 B 그래프용 데이터

현재 LLS(2열 덤프)로 한 단계 B 검증의 데이터와 그림입니다. 보고서 작성 지시는 `REPORT_BRIEF.md`에 있습니다.

| 위치 | 내용 | 만드는 방법 |
|---|---|---|
| `notes/*.csv` | 체크리스트 기록(`docs/checklist.html`의 사본)에서 뽑은 수치와 모델 값. 기록의 알려진 오류는 바로잡음 | `python3 lambda_patch/model/stageb_notes_to_csv.py` |
| `figures/*.png` | 보고서 그림 G1~G11 | `python3 lambda_patch/model/plot_stageb.py` (matplotlib) |
| (선택) 이 폴더의 `h0_*.csv`, `h1_*.csv` | MATLAB으로 다시 돌린 추가 곡선 (아래). 보고서 그림에는 필요 없고 지금은 없음 | MATLAB `export_stageb_csv` |

`notes/` 파일: `h0_frag.csv`(조건별 기존 Z·Λ̂ 오검출률, κ, Λ̂ 탭 꼬리), `h0_pkt.csv`(조건별 소프트 꼬리 비율·수락·쌍 상관), `tail_5k_vs_30k.csv`, `h1_levels.csv`(L90·L99와 모델), `h1_awgn_pd.csv`(AWGN 세기별 검출 수·Pd·95% 구간, 기존 Z·Λ̂ 프래그먼트와 Λ̂ 패킷 소프트·6-of-8, B-H1 메모 "G1 데이터"), `h1_pkt_pd_awgn_transition.csv`, `h1_mp_pd.csv`·`h1_mp_peak_offset.csv`·`h1_mp_firstpath.csv`(주거 LOS), `b0_peak_trajectories.csv`, `bmd_pass_count.csv`, `b1r_tail_summary.csv`·`b1r_tail_flagged.csv`·`b1r_packet_check.csv`, `consistency.csv`, `model_z_pfa.csv`. 각 행의 `source` 열에 출처가 있습니다.

## (선택) MATLAB 추가 곡선 (`export_stageb_csv`)

보고서 그림 G1~G11은 `notes/`만으로 그려집니다. 아래는 RIF 1·2·4개의 세기별 Pd, Λ̂ 탭 꼬리 곡선, 패킷 꼬리 곡선처럼 메모에 없는 추가 곡선이 필요할 때 씁니다.

### 만드는 방법

MATLAB에서 `lambda_patch`를 경로에 두고 실행합니다(덤프 루트 `D:/rif`, 출력은 이 폴더).

```matlab
export_stageb_csv                 % = export_stageb_csv('D:/rif', '<이 폴더>')
```

- 덤프 세트: `awgn` = `h0_awgn` / `h1_awgn`(AWGN, 0.25 ppm), `mp` = `h0_mp` / `h1_mp`(`ChanType=RESIDENTIAL LoSPath=1`, 802.15.4a 주거 LOS), `cfo0` = `h0_cfo0`(잔여 CFO 0). 없는 폴더는 건너뜁니다.
- 옵션은 `TEST_PLAN.md` 3.2절과 같습니다: Λ̂ 패킷은 `'Chip', 'fixed', 'Profile', 'fira-medium', 'Drift', 0`, 위치 정확도는 `'PeakTap', 125`(`'Tol'` AWGN 3, 멀티패스 2).
- H0는 `h0_stats_cache.mat`를 써서 빠르고, H1은 8번 모든 덤프를 다시 읽어 가장 오래 걸립니다.

### 파일

세기는 dBm(음수)입니다. 판정기: Z = 기존 Z 판정기, L = CIR 추정 Λ̂(프래그먼트).

| 파일 | 열 | 쓰는 그래프 |
|---|---|---|
| `h1_frag.csv` | set, len, dBm, n, Z_k, Z_Pd, Z_lo, Z_hi, Z_PosAcc, L_k, L_Pd, L_lo, L_hi, L_PosAcc (lo·hi: Clopper-Pearson 95%) | 감도 곡선(프래그먼트), 멀티패스 |
| `h1_pkt.csv` | set, rule, len, dBm, n, k, Pd, lo, hi, PosAcc. rule: `soft8`(Λ̂ 소프트, RIF 8개, 10⁻⁶), `kofn6`(6-of-8), `single1`·`soft2`·`soft4`(RIF 1·2·4개 구성) | 감도 곡선(패킷), RIF 구성, 멀티패스 |
| `h1_levels.csv` | set, rule, len, L90, L99, L90loc. rule은 `fragZ`, `fragL` 또는 위의 패킷 규칙. NaN은 스윕 안에서 닿지 않음 | 감도 요약 |
| `h0_frag.csv` | set, len, dBm, n, Z_k, Z_pfa, Z_upper, D_med_dB, kappa_med, L_k, L_pfa, L_upper, L_tail_t10, L_tail_t14, L_tail_t18 | 프래그먼트 오검출률 대 세기, κ |
| `h0_frag_tail.csv` | set, len, dBm, t, count, nTap, P_meas, P_chi2_2 (신호 창 탭에서 Λ̂ ≥ t, t = 0:0.25:40, 이론 e^(−t/2)) | Λ̂ 탭 꼬리 곡선 |
| `h0_pkt.csv` | set, K, len, dBm, packets, accepts, rhoPair, pairRatio, t1~t3, ratio1~ratio3 (P(χ²₂K ≥ t)가 1e-2 / 1e-3 / 1e-4 아래로 처음 내려가는 격자점의 꼬리 비율, 기대 개수가 5 미만이면 NaN) | 패킷 오검출 근거(꼬리 비율) |
| `h0_pkt_tail.csv` | set, K, len, dBm, t, count, paths, P_meas, P_chi2_2K (이동 0 소프트 합 경로에서 S ≥ t, 정수 격자 0:ceil(T)) | 패킷 꼬리 곡선 |

### 읽을 때 주의

- 패킷 수락은 모든 조건에서 0이라, 10⁻⁶ 자체는 셀 수 없습니다. 판정 근거는 꼬리 비율이고, 조건마다 5,000개(625패킷)에서는 실제 흔들림이 Poisson 오차의 1.3~1.5배입니다(`TEST_PLAN.md` 3.3절 "단계 B 결론").
- RIF 4개 −40 dBm처럼 경로 꼬리 비율이 깊어질수록 오르는 조건은, 경로 수가 아니라 패킷 수 ÷ 합집합 한계로 판단합니다("B-1R 결과").
- 위치 정확도는 첫 경로(탭 125) 기준입니다. 주거 LOS에서 0.42~0.46인 것은 채널 특성입니다("B-H1-3 결과").
