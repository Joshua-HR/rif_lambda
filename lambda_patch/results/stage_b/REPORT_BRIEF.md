# 단계 B 결과 보고서 작성 지시서

현재 LLS(C++ 링크 레벨 시뮬레이터, 2열 CIR 덤프)에서 수행한 RIF 판정기 검증(단계 B, 2026-11-16 완료)의 결과와 분석을 보고서로 정리하기 위한 지시서입니다. 이 저장소를 읽을 수 있는 다른 에이전트에게 다음처럼 전달합니다.

> `lambda_patch/results/stage_b/REPORT_BRIEF.md`를 읽고 지시대로 보고서를 만들어줘. 독자는 [예: 팀 내부 엔지니어 / 상위 보고], 결과물은 [예: 새 브랜치의 `lambda_patch/results/stage_b/REPORT.md`]야.

대괄호 두 곳은 맡기는 사람이 정합니다. 팀 내부용이면 5장 횡단 분석을 자세히, 상위 보고용이면 1장 요약과 6장 위험 위주로 짧게 씁니다.

## 1. 할 일

이 폴더의 데이터와 그림, 그리고 아래 원자료를 근거로 단계 B의 결과와 분석을 보고서로 정리합니다. 새 측정이나 새 계산으로 결론을 바꾸지 말고, 저장소에 기록된 결과를 독자가 따라갈 수 있게 정리하는 것이 목적입니다.

## 2. 원자료 (이것만 근거로 쓸 것)

| 우선 | 파일 | 볼 곳 |
|---|---|---|
| 1 | `lambda_patch/DECISIONS.md` | "결정 1: LLS 검증 (단계 B)"의 시험별 결과 표와 해석 항목, 결정 5 "피크 이동 후보"(LLS와 silicon의 구조 차이), 결정 6 "LLS 기준선"(B-F·B-S), 결정 3(지원 구성) |
| 1 | `lambda_patch/TEST_PLAN.md` | 1.2절(LLS 공통 인자, 멀티패스 채널), 3.1절(덤프 세트), 3.2절(MATLAB 명령), 3.3절(기대값과 결과: "B-0 결과", "B-MD 결과", "B-H0-3 결과", "B-H1-3 결과", "B-F·B-S 결과", "구성별 기대 감도", "B-1R 결과", "단계 B 결론", "단계 B 완료"), 8장, 9장(감시 항목, 첫 경로 탐색) |
| 2 | `lambda_patch/LAMBDA_PATCH_DESIGN.md` | 모델 근거와 모델 대비 LLS 비교: 2.1, 6.x, 7.4, 9.x절 |
| 2 | 이 폴더 `notes/*.csv`, `figures/*.png` | 체크리스트 기록에서 뽑은 수치와 그림(4장). 문서와 다르면 문서를 따름 |
| 3 | `docs/checklist.html`의 `window.RIF_CHECKLIST_SNAPSHOT` | 항목별 원자료 메모(MATLAB 출력, 2026-11-16 기준). 6장의 오류가 남아 있으니 `notes/*.csv`를 먼저 씀 |
| 3 | `lambda_patch/*.m` | 지표 정의를 확인할 때만 (`h0_pfa_analyze`, `h1_pd_analyze`, `rif_packet`, `rif_cir_analyze`, `check_packet_assumptions`) |

## 3. 보고서 구성

1. **요약 (반 쪽):** 결론 한 문장, 핵심 수치 5~7개, 남은 위험
2. **판정 기준과 대상**
   - 패킷 오수락 10⁻⁶ (FiRa medium, 측정 1회당)
   - 고정 칩 CIR 추정 Λ̂ 소프트 패킷 판정, 문턱은 목표/1.5
   - RIF 8개 = 프래그먼트 8개
   - 프래그먼트 1e-3과 기존 Z 판정기는 진단용
3. **LLS 설정**
   - HPRF, 잔여 오프셋 0.25 ppm, 2열 덤프
   - 채널: AWGN, `ChanType=RESIDENTIAL LoSPath=1`(802.15.4a 주거 LOS, τrms 중앙값 15.2 ns)
   - LLS 구조상 특징: 설정한 오프셋으로 RIF 시작점을 정수 샘플로 옮겨 피크가 1탭 톱니가 됨 → `'Drift'` 0
   - 피크 탭은 125
   - silicon의 오프셋 방식은 미정
4. **시험별 결과 표**
   - 대상 시험: B-0, B-CH, B-H0-4, B-H0, B-H1, B-F, B-S, B-MD, B-1R, B-R(B-H0-3, B-H1-2, B-H1-3, B-H0-5 생략)
   - 열: 목적, 덤프와 표본 수, 핵심 수치, 기준, 판정
5. **횡단 분석**
   - 소수 탭 톱니에 의한 소프트 합 손해 (L90 약 0.5 dB, L99 약 1 dB, RIF 개수별 차이)
   - 5,000개 조건의 꼬리 비율 흔들림이 Poisson의 1.3~1.5배인 이유와 30,000개 연장 사례
   - 경로 수가 아니라 패킷 수 ÷ 합집합 한계로 판단하는 이유 (B-1R)
   - 기존 Z 판정기에 패킷 보장이 없는 근거: B-MD의 같은 채널 몰림, B-H0-3의 잔여 CFO 의존
   - 주거 LOS 멀티패스: 감도 손해, 위치 정확도 0.42~0.46과 그 원인
   - 모델 대 LLS 비교: 모델이 비관적인 정도, 모델 채널(고정 LOS)과 LLS 채널의 차이
6. **감시 항목과 남은 일**
   - 감시 항목 2개와 확정 방법
   - B-H0-5 생략 근거
   - 첫 경로 탐색 과제
   - 단계 C·D(5열 덤프, 원시 샘플 Λ)
   - silicon 오프셋 방식
7. **부록:** 시험별 MATLAB 명령 위치(TEST_PLAN 절 번호), 데이터 파일과 그림 목록

## 4. 그래프

그림은 `figures/`에 있고, 다시 그리려면 5장의 방법을 씁니다. 본문에는 1·2순위를 반드시 넣고, 캡션에 근거 수치의 출처를 적습니다. 그림 라벨은 영어입니다(한글 글꼴을 가정하지 않음). 보고서 본문은 한국어로 쓰고, 필요하면 같은 데이터로 다시 그려도 됩니다.

| # | 그림 파일 | 내용 | 데이터 | 상태 |
|---|---|---|---|---|
| G1 | `g1_packet_pd_awgn_transition.png` (전체판: `g1_sensitivity_curves.png`) | 감도 곡선: Pd 대 입력 세기 | `notes/h1_pkt_pd_awgn_transition.csv`, 전체판은 `h1_frag.csv`·`h1_pkt.csv` | **일부.** 지금은 패킷 곡선의 전이 구간만 있음. 전체판(Λ̂ 패킷, 6-of-8, 기존 Z, Λ̂ 프래그먼트와 95% 구간)은 MATLAB 내보내기 뒤 생김 |
| G2 | `g2_sensitivity_by_configuration.png` | 감도 요약: L90·L99 대 RIF 길이, RIF 1·2·4·8개, 측정 대 모델 | `notes/h1_levels.csv` | 완료 |
| G3 | `g3_fragment_false_alarm.png` | 프래그먼트 오검출률 대 입력 세기: 기존 Z(측정·모델) 대 Λ̂ | `notes/h0_frag.csv`, `notes/model_z_pfa.csv` | 완료 |
| G4 | `g4_packet_tail_ratio.png` | 패킷 오검출 근거: 조건별 소프트 꼬리 비율, 0.7~1.3 띠, 수락 0 | `notes/h0_pkt.csv` | 완료 (꼬리 곡선 P(S ≥ t)는 `h0_pkt_tail.csv`가 생기면 추가 가능) |
| G5 | `g5_peak_sawtooth.png` | LLS 피크 톱니 (B-0 예시 패킷) | `notes/b0_peak_trajectories.csv` | 완료 |
| G6 | `g6_softsum_margin.png` | 소프트 합 손해 대 RIF 개수 (모델 − 측정) | `notes/h1_levels.csv` | 완료 |
| G7 | `g7_z_detector_pass_count.png` | 기존 Z의 패킷당 통과 개수 대 이항 (B-MD) | `notes/bmd_pass_count.csv` | 완료 |
| G8 | `g8_multipath.png` | 주거 LOS: 감도, 위치 정확도, 피크 오프셋 (B-H1-3) | `notes/h1_mp_pd.csv`, `h1_mp_peak_offset.csv`, `h1_mp_firstpath.csv` | 완료 |
| G9 | `g9_residual_cfo.png` | 잔여 CFO 0.25 ppm 대 0: κ와 기존 Z 오검출 | `notes/h0_frag.csv` | 완료 (부록) |
| G10 | `g10_rif4_watch_item.png` | B-1R 감시 항목: 경로·패킷 초과 ÷ 합집합 한계 대 깊이 | `notes/b1r_packet_check.csv` | 완료 (부록) |
| G11 | `g11_5k_vs_30k.png` | 5,000개 대 30,000개 꼬리 비율 | `notes/tail_5k_vs_30k.csv` | 완료 (부록) |
| 표 | — | 일관성 검사(B-F·B-S), RIF 4개 꼬리 요약 | `notes/consistency.csv`, `b1r_tail_summary.csv`, `b1r_tail_flagged.csv` | 표로 넣음 |

## 5. 데이터와 그림을 다시 만드는 방법

| 단계 | 명령 | 만드는 것 |
|---|---|---|
| 체크리스트 기록 → CSV | `python3 lambda_patch/model/stageb_notes_to_csv.py` | `notes/*.csv` (6장의 오류를 코드에서 바로잡음) |
| MATLAB 전체 곡선 (사람이 실행) | MATLAB: `export_stageb_csv` (덤프 `D:/rif`) | `h1_frag.csv`, `h1_pkt.csv`, `h1_levels.csv`, `h0_frag.csv`, `h0_frag_tail.csv`, `h0_pkt.csv`, `h0_pkt_tail.csv` (`README.md`) |
| 그림 | `python3 lambda_patch/model/plot_stageb.py` (matplotlib 필요) | `figures/*.png`. MATLAB CSV가 있으면 G1 전체판을 그림 |

- 다른 에이전트는 MATLAB과 덤프에 접근할 수 없습니다. 필요한 데이터가 저장소에 없으면 만들어 내지 말고 "데이터 필요"로 표시하고 무엇이 필요한지 적어 주세요.
- G1 전체판, 꼬리 곡선, 다른 조건이 필요하면 맡긴 사람이 MATLAB에서 `export_stageb_csv`를 돌려 CSV를 이 폴더에 커밋합니다.
- 체크리스트에 새 결과가 적히면 Pages 사본(`docs/checklist.html`)을 먼저 갱신한 뒤 `stageb_notes_to_csv.py`를 다시 돌립니다.

## 6. 지켜야 할 것

- **출처:** 수치는 2장의 원자료에 있는 것만 쓰고, 각 수치 옆에 출처(파일과 절, 항목 ID 또는 CSV 파일)를 붙입니다. 새로 계산하거나 추정한 값은 "추정"으로 따로 표시합니다.
- **모델과 측정:** "모델 기대값"과 "LLS 측정값"을 항상 구분합니다.
- **원자료 메모의 오류:** `docs/checklist.html` 사본에는 문서에서 바로잡은 오류가 남아 있습니다. 문서와 `notes/*.csv`를 따르고, 메모를 직접 읽을 때도 다음 값은 쓰지 않습니다.
  - B-H0-4 메모의 "80.0"은 −82 dBm입니다. `h0_mp`에는 −40·−82 dBm만 있습니다.
  - B-1R 메모 (1)의 세 줄은 같은 값으로 섞인 출력입니다. "꼬리 비율 재출력 결과"를 씁니다.
  - B-1R "패킷 초과 결과"의 "기대" 열은 약 50배 작게 잘못 계산된 값입니다. 맞는 값은 `notes/b1r_packet_check.csv`와 TEST_PLAN 3.3 "B-1R 결과"에 있습니다.
  - B-H0 메모 (2)는 1e3 배율로 출력되어 꼬리 비율 해상도가 0.1입니다. 64심볼 30,000개 표의 "60 120"은 64 / −120입니다.
  - B-H1 메모 (2)의 위치 정확도 열은 예전 기준(탭 127 ± 2)이라 쓰지 않습니다. 위치 정확도는 탭 125 기준 1.00입니다.
- **용어:** 문서와 같게 씁니다. "원시 샘플 Λ", "CIR 추정 Λ̂", "기존 Z 판정기", "다상 성분", "위상 간 바닥비".
- **날짜:** 절대 날짜(예: 2026-11-16)로 씁니다.
- **고치지 말 것:** `RIF_CIR_Validity_Detector.md`, `.m` 파일, 기존 문서(DECISIONS, TEST_PLAN, 설계 문서)는 고치지 않습니다. 보고서와 새 그림만 추가합니다.
- **브랜치:** master에 직접 커밋하지 않습니다. 새 브랜치에서 작업하고, 끝나면 브랜치 이름을 알려 줍니다.
- **모르는 것:** 확인이 안 되는 내용은 추측하지 말고 "문서에 없음"으로 남깁니다.

## 7. 끝나면 알려 줄 것

- 만든 파일과 브랜치 이름
- 보고서에 쓴 핵심 수치와 그 출처 목록 (검토용)
- "데이터 필요"로 남긴 그래프나 표, 그리고 무엇이 있으면 채울 수 있는지
- 원자료끼리 맞지 않아 판단을 보류한 곳
