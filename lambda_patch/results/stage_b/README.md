# 단계 B 그래프용 CSV

현재 LLS(2열 덤프)로 한 단계 B 검증의 곡선 데이터입니다. 체크리스트 메모와 `TEST_PLAN.md` 3.3절에는 요약값만 있으므로, 보고서 그래프는 이 파일들로 그립니다.

## 만드는 방법

MATLAB에서 `lambda_patch`를 경로에 두고 실행합니다(덤프 루트 `D:/rif`, 출력은 이 폴더).

```matlab
export_stageb_csv                 % = export_stageb_csv('D:/rif', '<이 폴더>')
```

- 덤프 세트: `awgn` = `h0_awgn` / `h1_awgn`(AWGN, 0.25 ppm), `mp` = `h0_mp` / `h1_mp`(`ChanType=RESIDENTIAL LoSPath=1`, 802.15.4a 주거 LOS), `cfo0` = `h0_cfo0`(잔여 CFO 0). 없는 폴더는 건너뜁니다.
- 옵션은 `TEST_PLAN.md` 3.2절과 같습니다: Λ̂ 패킷은 `'Chip', 'fixed', 'Profile', 'fira-medium', 'Drift', 0`, 위치 정확도는 `'PeakTap', 125`(`'Tol'` AWGN 3, 멀티패스 2).
- H0는 `h0_stats_cache.mat`를 써서 빠르고, H1은 8번 모든 덤프를 다시 읽어 가장 오래 걸립니다.

## 파일

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

## 읽을 때 주의

- 패킷 수락은 모든 조건에서 0이라, 10⁻⁶ 자체는 셀 수 없습니다. 판정 근거는 꼬리 비율이고, 조건마다 5,000개(625패킷)에서는 실제 흔들림이 Poisson 오차의 1.3~1.5배입니다(`TEST_PLAN.md` 3.3절 "단계 B 결론").
- RIF 4개 −40 dBm처럼 경로 꼬리 비율이 깊어질수록 오르는 조건은, 경로 수가 아니라 패킷 수 ÷ 합집합 한계로 판단합니다("B-1R 결과").
- 위치 정확도는 첫 경로(탭 125) 기준입니다. 주거 LOS에서 0.42~0.46인 것은 채널 특성입니다("B-H1-3 결과").
