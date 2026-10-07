# RIF Λ 판정기 검증 계획

- 작성일: 2026-10-07
- 대상: `lambda_patch/`의 MATLAB 분석 패치(`rif_cir_analyze`, `h0_pfa_analyze`, `h1_pd_analyze`)와 LLS(C++) RIF CIR 덤프
- 관련 문서: 설계 `LAMBDA_PATCH_DESIGN.md`, 해설 페이지 `../docs/index.html`
- 이 문서의 기대값은 모델(`model/`)과 합성 덤프에서 나온 값입니다. MATLAB과 LLS 실측으로 확인하기 전까지는 기대값으로만 보세요.

---

## 0. 개요

| 단계 | 목적 | C++ 수정 | 데이터 |
|---|---|---|---|
| **A** | MATLAB 도구 점검. **5열과 2열 합성 덤프를 MATLAB으로 각각 만들어 확인** | 없음 | 합성 5열 / 2열 |
| **B** | 현재 LLS로 기준선 재현, MD 방식 실패 확인, **고정 칩 대안(Λ̂) 결론** | 없음 | LLS 2열 |
| **C** | Q·Π 덤프 구현 확인 (bring-up) | 있음 | LLS 5열 소량 |
| **D** | 원시 샘플 Λ 최종 검증 | 있음 | LLS 5열 |
| E (선택) | 거리 ≠ 0, 고정소수점, 다중 안테나, BPRF | 상황별 | LLS 5열 |

단계 A의 2열 세트는 단계 B(현재 LLS)와 형식이 같습니다. 그래서 B에서 쓸 명령과 결과 표를 LLS를 돌리기 전에 미리 연습할 수 있습니다.

## 1. 공통 준비

### 1.1 MATLAB

```matlab
% 원본 스크립트가 있는 저장소 최상위가 아닌 폴더에서 실행 (현재 폴더의 함수가 path보다 우선)
addpath('<repo>/lambda_patch', '-begin');
addpath('<repo>/lambda_patch/test');
which rif_cir_analyze      % <repo>\lambda_patch\rif_cir_analyze.m 이어야 함
```

### 1.2 LLS 공통 설정 (MD 문서 11장)

```
Mode=HPRF Rate=R6M MMSEn=1 MMSOnlyEn=1 NumRSF=0 NumRIF=4 MMRSSeq=33 MMRSGap=33 NumRx=1
AdcBit=6 TargetRms=31 CfoEn=1 ScoEn=1 CfoEstEn=0 ScoTrkEn=0 ScoLinkedToCFO=0
SemiSyncOnlyMode=1 MMSPwCtrlOn=1  FreqOffset=0.25 ClkOffset=0.25  ChanType=AWGN
RIFStsLen=0..3 (32/64/128/256심볼), SemiSyncAcclen = 심볼 수
```

- **H0(틀린 키):** Rx STS 시드만 Tx와 다르게 둡니다.
- **H1(맞는 키):** Tx와 Rx 시드를 같게 둡니다.
- **job 독립성:** job마다 STS 시드와 잡음 시드를 모두 다르게 해야 합니다. 같으면 `h0_pfa_analyze`가 "identical Z values" 경고를 냅니다.

### 1.3 덤프 폴더 규칙 (합성·LLS 공통)

```
<root>/c<심볼수>_m<세기 절댓값, 정수>_j<job>/bin/<SignalPower>_RifCir_AccNum_<N>_Frame<f>_Samp<s>.txt
```

- **덤프 켜기 (LLS):** `semisync_only()`의 `if (0 && iSemiSyncNextState_ == RX_STS)`를 `1`로 바꿉니다.
- **실행 위치:** job마다 실행 폴더를 따로 둡니다. `uwb/bin`은 피하세요(덤프가 git 추적 대상으로 잡힘).
- **root 분리:** H0와 H1, 채널, CFO 조건은 root를 따로 둡니다.
- **프레임 수:** 프레임당 8 프래그먼트이므로 프레임 수 = 원하는 프래그먼트 수 ÷ 8입니다.

### 1.4 표본 수

| 용도 | 조건당 프래그먼트 | 근거 |
|---|---|---|
| H0 선별 | 5,000 | 탭 꼬리 비율 표본 28.5만 개 |
| H0 최종 판정 | 30,000 | 단측 95% 상한 ≤ 1e-3. n = 5,000이면 참 Pfa 5e-4도 합격 확률 29% |
| H1 | 세기당 1,000 (L99는 문턱 근처 2,000) | Pd 0.9 근처 신뢰구간 ±2% |
| 합성 리허설 | 세기당 500~1,000 | 파이프라인 확인용이라 통계보다 속도 우선 |

덤프 크기는 대략 2열 5~7 KB, 5열 13~15 KB입니다. LLS에서 한 조건을 먼저 돌려 크기와 프레임당 시간을 재세요.

## 2. 단계 A: MATLAB 도구 점검 (LLS 불필요)

### A-1. 5열 합성 세트

```matlab
gen_synthetic_dumps('D:/rif_test5');            % h0, h0_mp20, h1, h0_2col, scale_test, mixed_test (7,620개)
res5 = run_synthetic_test('D:/rif_test5');      % 5열 자동 인식 → 9단계, 43개 항목
```

### A-2. 2열 합성 세트 (현재 LLS 형식, 같은 시드)

```matlab
gen_synthetic_dumps('D:/rif_test2', 'Cols', 2); % h0, h0_mp20, h1만 (5열 전용 폴더 3개는 건너뜀)
res2 = run_synthetic_test('D:/rif_test2');      % 2열 자동 인식 → 5단계, 29개 항목
```

| 2열 단계 | 확인 내용 |
|---|---|
| 1 | 폴더 하나: Q·Π 없음(`hasMom` false), raw Λ 미계산, MD 오검출, κ, `'Moments','cir'`로 Λ̂ 계산과 위상당 추정 탭 수(24~25) |
| 2 | H0 전체: raw Λ 표 없음, 세기별 MD Pfa와 κ |
| 3 | H0 전체 `'cir'`: Λ̂ 표, 최대 Pfa, 탭 꼬리 비율 |
| 4 | 멀티패스: D 중앙값 < 10 dB, MD Pfa, Λ̂ Pfa |
| 5 | H1: MD L90, raw Λ의 Pd가 0, Λ̂ L90 − MD L90, 위치 정확도, −110 dBm 신뢰구간 |

### A-3. 2열·5열 교차 확인 (정확히 같아야 함)

같은 프리셋과 시드로 만들었으므로 두 세트의 CIR은 비트 단위로 같습니다. 그래서 CIR만 쓰는 결과는 근사가 아니라 정확히 일치해야 합니다.

```matlab
R5  = h0_pfa_analyze('D:/rif_test5/h0', 'Suggest', false, 'Plot', false);
R2  = h0_pfa_analyze('D:/rif_test2/h0', 'Suggest', false, 'Plot', false);
isequal(R5.E.k, R2.E.k)                                  % MD 오검출 수        → true
isequal([R5.cond.kappa], [R2.cond.kappa])                % kappa               → true
R5c = h0_pfa_analyze('D:/rif_test5/h0', 'Suggest', false, 'Moments', 'cir', 'Plot', false);
R2c = h0_pfa_analyze('D:/rif_test2/h0', 'Suggest', false, 'Moments', 'cir', 'Plot', false);
isequal(R5c.EL.k, R2c.EL.k)                              % Λ̂ 오검출 수         → true
R15 = h1_pd_analyze('D:/rif_test5/h1', 'Moments', 'cir', 'Plot', false);
R12 = h1_pd_analyze('D:/rif_test2/h1', 'Moments', 'cir', 'Plot', false);
isequal(R15.table, R12.table)                            % H1 Pd 표 (MD, Λ̂)   → true
```

하나라도 `false`면 덤프 읽기나 형식 처리에 문제가 있는 것입니다.

결과를 그림으로 나란히 보려면 JSON으로 내보내 웹 페이지 9장 "1단계 결과: 5열 덤프 대 2열 덤프" 그림에서 불러옵니다.

```matlab
export_web_results('D:/rif_test5', 'D:/rif_test2', 'D:/rif_results.json', 'Res5', res5, 'Res2', res2);
```

- 그림 위의 "MATLAB 결과 불러오기 (JSON)" 버튼으로 파일을 고릅니다. 파일은 브라우저 안에서만 읽습니다.
- 교차 확인 타일 4개(MD 오검출 수, κ, Λ̂ 오검출 수, H1 검출 수)가 모두 "같음"이어야 합니다.
- `run_synthetic_test` 타일에는 43/43 · 29/29가 보여야 합니다.
- 미리 만든 `synthetic_dumps/`와 `synthetic_dumps_2col/`로 내보내면 "Python 미러와 비교" 타일이 함께 나옵니다. 경계값 몇 건 외에는 "같음"이어야 합니다.
- GitHub Pages 판은 `docs/results/matlab_results.json`이 있으면 그것을 기본으로 보여 줍니다.

### A-4. 단계 B 리허설 (선택)

단계 B와 같은 폴더 구성을 합성 2열로 만들어, LLS를 돌리기 전에 명령·폴더·결과 표를 한 번 돌려 봅니다.

```matlab
root = 'D:/rif_rehearsal';  sd = 100;
for p = [120 100 95 90 88 86 84 82 80 78 76 72 60 40]            % B-H0-1
    for j = 1:2
        sd = sd + 1;
        gen_synthetic_dumps(fullfile(root, 'h0_awgn', sprintf('c32_m%d_j%d', p, j), 'bin'), 'Preset', 'none', ...
                            'Cols', 2, 'Power', -p, 'Key', 'wrong', 'Count', 500, 'Seed', sd);
    end
end
for p = 100:120                                                  % B-H1-1
    sd = sd + 1;
    gen_synthetic_dumps(fullfile(root, 'h1_awgn', sprintf('c32_m%d_j1', p), 'bin'), 'Preset', 'none', ...
                        'Cols', 2, 'Power', -p, 'Key', 'right', 'Count', 500, 'Seed', sd);
end
```

- 다른 조건은 `'N', 64|128|256`, `'Ppm', 0`, `'Trms', 20, 'KdB', -6`으로 만듭니다.
- 이 결과는 단계 B의 "모델 기대값"과 같은 출처이므로, 나중에 LLS 결과와 나란히 놓고 비교하면 모델과 LLS의 차이가 바로 보입니다.

### 단계 A 기대값과 합격 기준

| 항목 | A-1 (5열) | A-2 (2열) |
|---|---|---|
| 요약 | `43 / 43 checks OK` | `29 / 29 checks OK` |
| MD Pfa (−120/−88/−82/−40 dBm) | 약 0.1% / 1.0% / 5.2% / 0.2% | 같음 |
| κ 중앙값 | 0.38 / 0.57 / 0.82 / 0.98 | 같음 |
| Λ (raw) Pfa | 0~0.1% | 해당 없음 |
| Λ̂ (cir) Pfa | 0~0.2% | 같음 |
| H1 L90 | MD −111.2, Λ −112.3, Λ̂ −109.7 dBm | MD −111.2, Λ̂ −109.7 dBm |

이 숫자는 미리 만든 `synthetic_dumps/`의 값입니다. MATLAB에서 새로 만든 세트는 난수가 달라 조금 다르며, 기대 범위는 그 흔들림을 감안했습니다. A-3의 교차 확인은 모두 `true`여야 합니다.

## 3. 단계 B: 현재 LLS (2열)

### 3.1 덤프 세트

| ID | 가설 | 길이 | 세기 [dBm] | 채널 / CFO | 조건당 n | 우선 |
|---|---|---|---|---|---|---|
| B-H0-1 | H0 | 32 | −120, −100, −95, −90, −88, −86, −84, −82, −80, −78, −76, −72, −60, −40 | AWGN / 0.25 ppm | 5,000 | **필수** |
| B-H0-2 | H0 | 64, 128, 256 | −120, −82, −40 | AWGN / 0.25 ppm | 5,000 | **필수** |
| B-H0-3 | H0 | 32, 256 | −82, −40 | AWGN / **FreqOffset=0, ClkOffset=0** | 5,000 | 권장 |
| B-H0-4 | H0 | 32 | −82, −40 | **멀티패스** (LLS 실내 채널 1~2종) | 5,000 | 권장 |
| B-H0-5 | H0 | 32 | Tx off (불가하면 −150) | AWGN | 5,000 | 권장 |
| B-H1-1 | H1 | 32, 256 | −120 … −100, 1 dB | AWGN | 1,000 | **필수** |
| B-H1-2 | H1 | 64, 128 | −120 … −100, 1 dB | AWGN | 1,000 | 권장 |
| B-H1-3 | H1 | 32 | −120 … −98, 2 dB | B-H0-4의 멀티패스 | 1,000 | 권장 |

### 3.2 실행 (A-2·A-4와 같은 명령)

```matlab
R  = rif_cir_analyze('D:/rif/h0_awgn/c32_m82_j1/bin', 'Plot', 'both', 'Index', 1);     % 표본 점검
R0  = h0_pfa_analyze('D:/rif/h0_awgn', 'Suggest', false, 'Save', 'D:/rif/report/B_h0_md');
R0c = h0_pfa_analyze('D:/rif/h0_awgn', 'Suggest', false, 'Moments', 'cir', 'Save', 'D:/rif/report/B_h0_cir');
R1  = h1_pd_analyze('D:/rif/h1_awgn');                     % MD Pd·CI·위치 정확도·L90
R1c = h1_pd_analyze('D:/rif/h1_awgn', 'Moments', 'cir');   % Λ̂
saveas(gcf, 'D:/rif/report/B_h1_cir.png');
% Λ̂ 전제 확인 (대표 폴더): 추정 범위를 신호 창 앞쪽으로만 바꿔도 비슷해야 함
a = rif_cir_analyze(fld, 'Quiet', true, 'Plot', 'none', 'Moments', 'cir');
b = rif_cir_analyze(fld, 'Quiet', true, 'Plot', 'none', 'Moments', 'cir', 'CirNoise', [0 118]);
```

### 3.3 기대값과 판정 (모델, 32심볼; 모델은 MD 실측보다 약 0.5 dB 비관적)

| 확인 항목 | 볼 곳 | 기대값 | 판정 |
|---|---|---|---|
| MD 7.1 표 재현 | `R0`, −120 / −40 dBm | 약 1e-3 / 약 2e-3 | MD 실측과 같은 수준 |
| **MD 중간 세기 실패** | `R0`, −88 ~ −76 dBm | Pfa 1e-2 ~ 5.6e-2, D 중앙값 < 10 dB | 모델 예측 재현 여부 |
| κ (비원형 원인) | `R0` kappa 열 | −120: 0.38 / −40 dBm: 32심볼 0.98, 64심볼 0.91, 128심볼 0.66, 256심볼 0.38 | 길이에 따른 감소 모양 일치 |
| CFO 0 | B-H0-3 | 256심볼도 κ ≈ 1, MD Pfa ≈ 2.7e-3 | |
| 멀티패스 | B-H0-4 | D < 10 dB 다수, MD Pfa 3~5e-2 | |
| **Λ̂ 오검출** | `R0c` Λ 표 | 모든 조건 ≤ 1e-3 | 점추정 ≤ 1e-3 (n = 5,000에서 k ≤ 5) |
| **Λ̂ 탭 꼬리 비율** | `R0c` ratio t = 10, 14 | 0.8 ~ 1.2 | **0.7 ~ 1.3 → 전제 성립** |
| Λ̂ 감도 | `R1c.L90` | MD −111.3 (실측 −111.9), Λ̂ −109.6 (LLS 약 −110.1 예상) | MD 대비 손해 ≤ 2.5 dB |
| 위치 정확도 | `R1c.cond(c).PosAcc` | AWGN ≥ 98% | AWGN ≥ 98% |
| 합성과의 비교 | A-4 결과 vs B 결과 | 표가 비슷한 모양 | 크게 다르면 모델과 LLS 차이(ADC, AGC 등) 기록 |

**단계 B 결론 기준:** Λ̂의 H0 Pfa가 모든 조건에서 1e-3 이하이고 탭 꼬리 비율이 0.7~1.3이면 "고정 칩에서 Λ̂ 사용 가능"으로 정리합니다.

## 4. 단계 C: 5열 덤프 bring-up (LLS 소량)

C++ 변경은 설계 문서 4.2절과 같습니다(`momQ += |z|²`, `momPi += z²`, 같은 z, 같은 AccLen−1 구간, 프래그먼트마다 초기화).

| ID | 조건 (조건당 200~500개) | 확인 | 기준 |
|---|---|---|---|
| C-1 | H0 −120 dBm | `median(R.qRatio)`, `median(R.rhoPk)` | 0.9~1.1, ≤ 0.1 |
| C-2 | H0 −40 dBm | `qRatio`, `rhoPk` | ≈ 1, ≈ 0.97 |
| C-3 | H1 −60 dBm | `validL`, `kpk` | 전부 1, 127 근처 |
| C-4 | C-1과 같은 시드로 2열/5열 | 1·2열(CIR) 비교 | 동일 (A-3과 같은 방식) |

고장 진단 (아래 외에, 같은 시드인데 C-4에서 CIR이 다르면 Q·Π 누적이 상관 계산 경로를 바꾼 것입니다)

| 증상 | 원인 후보 |
|---|---|
| `qRatio`가 1/64 또는 64 근처 (`rif_cir_analyze:scale` 경고) | C와 Q·Π의 라운딩 스케일 불일치 |
| `qRatio`가 1에서 10~30% 벗어남 | Q를 다른 샘플·구간에서 누적 |
| 열잡음에서 `rhoPk` ≈ 1 | Π를 `z·conj(z)`로 계산 |
| `rif_cir_analyze:mixed` 오류 | 한 폴더에 이전 2열 덤프가 남아 있음 |

## 5. 단계 D: 5열 최종 검증

| ID | 내용 | n |
|---|---|---|
| D-H0-핵심 | 32심볼 −120, −88, −82, −76, −40 dBm / 256심볼 −82, −40 dBm (CFO 0 포함) / 멀티패스 −40 dBm | **30,000** |
| D-H0-나머지 | 단계 B의 나머지 H0 | 5,000 |
| D-H1 | 단계 B의 H1 | 1,000 (문턱 근처 2,000) |

```matlab
R0  = h0_pfa_analyze('D:/rif5/h0_awgn', 'Suggest', false, 'Save', 'D:/rif/report/D_h0');   % raw Λ
R0h = h0_pfa_analyze('D:/rif5/h0_awgn', 'Suggest', false, 'PerClass', true);              % HW 근사
R0c = h0_pfa_analyze('D:/rif5/h0_awgn', 'Suggest', false, 'Moments', 'cir');              % 같은 덤프로 Λ̂
R1  = h1_pd_analyze('D:/rif5/h1_awgn');
R1c = h1_pd_analyze('D:/rif5/h1_awgn', 'Moments', 'cir');
```

| 합격 기준 | 기준 | 모델 기대값 |
|---|---|---|
| H0 탭 꼬리 비율 (t = 10, 14) | 0.7 ~ 1.3 | 0.83 ~ 1.25 |
| H0 핵심 조건 Pfa (n = 30,000) | Λ 표 `PASS` (단측 95% 상한 ≤ 1e-3) | 2.5e-4 ~ 8.5e-4 |
| H1 L90 (raw Λ) | MD 규칙보다 나쁘지 않을 것 | MD보다 약 0.6 dB 좋음 |
| H1 L90 (Λ̂, 같은 덤프) | raw 대비 손해 1.5 ~ 3.5 dB | 약 2.5 dB |
| 점수 방식 대비 (MD 기준 2) | 1 dB 이상 좋을 것 | +3.0 / +5.5 / +7.4 / +6.8 dB |
| PerClass | raw와 통계 오차 범위 안 | 차이 없음 |
| 위치 정확도 (AWGN) | ≥ 98% | 98 ~ 100% |

## 6. (선택) 단계 E

| 항목 | 방법 | 볼 것 |
|---|---|---|
| 거리 ≠ 0 | `ToaTestDist`. 분석은 `'Sig'`를 피크 주변으로, `'PeakTap'`을 실제 피크로 | 문턱은 \|W_s\| 공식대로 |
| 고정소수점 | `FXP_RX`가 정의되는 빌드에서 D-H0 핵심 조건 반복 | 탭 꼬리 비율 유지, 포화 없음 |
| 다중 안테나 | NumRx > 1, 기준 안테나 덤프 | 단일 안테나와 같은 H0 분포 |
| BPRF | `'Period', 16` | 같은 검사 |

## 7. 기록할 결과

| 결과 | 출처 |
|---|---|
| 단계 A 요약 (43/43, 29/29), A-3 교차 확인 결과 | `run_synthetic_test` 출력, `isequal` 결과, `export_web_results` JSON과 웹 9장 그림 |
| H0 MD 표 (Pfa, 상한, kappa, Dmed) | `h0_pfa_analyze` "current thresholds" 표, `h0_pfa.png` |
| H0 Λ 표 (raw / cir) | "Lambda detector" 표, `h0_lambda.png` |
| H1 표 (Pd [95% CI], 위치 정확도, L90/L99, L90loc) | `h1_pd_analyze` 출력과 그림 |
| bring-up 지표 (qRatio, rhoPk) | 단계 C의 `rif_cir_analyze` |
| 합성 대 LLS 비교 | A-4와 B의 같은 표를 나란히 |

결과 변수(`R0.cond`, `R0.EL`, `R1.cond`, `R1.L90` 등)는 `.mat`으로 저장해 두면 단계 간 비교가 쉽습니다.

## 8. 진행 체크리스트

1. [ ] A-1 5열 합성: 43/43
2. [ ] A-2 2열 합성: 29/29
3. [ ] A-3 교차 확인: 모두 `true`
4. [ ] (선택) A-4 단계 B 리허설
5. [ ] B 필수 세트(B-H0-1, B-H0-2, B-H1-1) → MD 실패 재현, κ, **Λ̂ 결론**
6. [ ] B 권장 세트(CFO 0, 멀티패스, Tx off, 나머지 길이)
7. [ ] C: C++ 5열 구현, bring-up 4항목
8. [ ] D: 5열 전체, 핵심 H0 30,000개
9. [ ] (선택) E

## 9. 보류 항목

- MMS 프래그먼트 소프트 결합(ΣΛ, 자유도 2K인 카이제곱 문턱)의 분석 도구 반영
- 첫 경로 탐색(ToA)과 위치 정확도의 연계
- RTL 비트폭, Q·Π 데시메이션
