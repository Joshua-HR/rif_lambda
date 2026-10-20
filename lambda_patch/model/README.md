# 통계·파형 모델 (Python 3, 표준 라이브러리만)

`LAMBDA_PATCH_DESIGN.md`와 `docs/index.html`에 인용한 수치를 다시 만드는 스크립트입니다.
numpy 없이 표준 라이브러리와 `multiprocessing`(6 프로세스)만 씁니다. 이 폴더에서 실행합니다.

| 파일 | 내용 | 재현하는 수치 | 시간(6코어) |
|---|---|---|---|
| `rifsim.py` | 통계 모델 기본: 빗 프로파일, CFO 공분산, MD 규칙 | 공용 모듈 | – |
| `rifsim3.py` | 멀티패스 채널(LOS + 지수 감쇠 확산) | 공용 모듈 | – |
| `rifsim4.py` | 통계 모델 + Q·Π + Λ (`run()`) | 공용 모듈 | – |
| `sweep_for_page.py` | 틀린 키 세기 스윕, MD 규칙 대 Λ (32심볼, 조건당 60,000) | 설계 문서 6.1 표, 웹 R1 그래프 | 약 3분 |
| `conds_for_page.py` | 9개 H0 조건 × 3 방식(MD / 정규화만 / Λ) | 설계 문서 6.2 표, 웹 R2 그래프 | 약 2분 |
| `l90_for_doc.py` | 길이별 L90: MD 규칙, MD 안전 단일 문턱, Λ | 설계 문서 6.4 표 | 약 1.5분 |
| `h1_pos_for_page.py` | H1 Pd와 95% 신뢰구간, 검출 위치 정확도 (AWGN, 세기당 10,002) | 설계 문서 6.5 표, 웹 R4 왼쪽 | 약 1분 |
| `h1_pos_mp_for_page.py` | H1 검출 위치 정확도, 멀티패스 τ10 / τ20 (세기당 4,002) | 설계 문서 6.5, 웹 R4 오른쪽 | 약 3분 |
| `cir_only_lambda.py` | 고정 칩 대안: CIR에서 Q·Π를 추정한 Λ̂ (H0 9조건 × 30,000, H1 세기당 6,000) `python3 cir_only_lambda.py [out.json]` | 설계 문서 8.3, 웹 8장 | 약 10분 |
| `wavesim.py` | 파형 수준 검증(실제 ±1 STS, 원시 샘플 Q·Π): `python3 wavesim.py 24000` | 설계 문서 6.3 표 | 약 15분 |
| `gen_dumps.py` | 합성 5열 덤프 생성: `python3 gen_dumps.py <root> [--N 64] [--ppm 0.25] [--drift D \| --sco S] [--frag-interval 1e-3]` (`<root> h1_drift`: 그 폴더만 추가). 이동량 기본은 잔여 CFO와 연동(0.2496탭), H1·이동 폴더 세기는 길이·CFO에 맞춰 옮김, 루트에 생성 조건 파일 `synthetic_set.txt` | `synthetic_dumps/` (`--drift 0.3`) | 약 6분(32심볼), 길이에 비례 |
| `verify_dumps.py` | MATLAB 패치 계산의 파이썬 판: `python3 verify_dumps.py <root>` | `run_synthetic_test` 기대값 | 수 초 |
| `make_2col.py` | 5열 세트에서 같은 CIR의 2열 세트 만들기: `python3 make_2col.py ../../synthetic_dumps ../../synthetic_dumps_2col` | `synthetic_dumps_2col/` | 수 초 |
| `web_results.py` | 5열·2열 세트 분석 결과 JSON (`export_web_results.m`과 같은 형식): `python3 web_results.py ../../synthetic_dumps ../../synthetic_dumps_2col out.json` | 웹 10장 "1단계 결과" 기본 데이터 | 약 10초 |
| `packet.py` | `rif_packet.m`과 `check_packet_assumptions.m`의 파이썬 판 (패킷 판정, 이동 기울기, AGC, 독립성) | `run_synthetic_test` 10단계 기대값 | 수 초 |
| `medium_level.py` | FiRa medium(패킷 10⁻⁶) 결합 규칙별 패킷 Pd와 H0 근거: `KF=8 python3 medium_level.py h1 4000 -122 -106 out.json`, `... h0 4000 out.json` | 설계 문서 9.3 첫 표 | 약 3분 + 1분 |
| `policy_k8.py` | 8개 RIF 정책 비교 (strict, soft, 일관성 하한): `python3 policy_k8.py 3000 -121 -104 out.json` | 설계 문서 9.3 정책 표 | 약 2분 |
| `single_len.py` | RIF 1개 설정, 길이별 L90 (1e-3, 10⁻⁶): `python3 single_len.py 3000 out.json` | 설계 문서 9.3 마지막 표 | 약 1분 |
| `shared_channel_md.py` | 패킷의 8개 프래그먼트가 멀티패스 채널 하나를 공유할 때(τrms 20 ns, 확산 6 dB 우세, 틀린 키 −42.5 / −82.5 dBm): MD 규칙의 채널별 프래그먼트 오검출 p, 8개 모두 통과의 패킷 오수락 E[p⁸], Λ̂ 소프트·AND 수락 수 (채널 1,998개, 9,990패킷) | 설계 문서 2.1절, 기준 전환 페이지 7장 | 미측정 |
| `shared_channel_tail.py` | 같은 조건에서 이동 0 소프트 합의 꼬리 비율 P(S ≥ t)/P(χ²₁₆ ≥ t), P = 1e-2 / 1e-3 / 1e-4 지점, 채널 공유와 프래그먼트별 채널 비교 (조건당 6,000패킷, 경로 34.2만 개) | 설계 문서 9.6·9.7절, 기준 전환 페이지 6장 | 미측정 |
| `length_l90.py` | RIF 길이별(32·64·128·256) 프래그먼트·패킷 규칙 10개의 L90, 다른 잔여 CFO에서 CFO 이득 공식 확인: `KF=8 python3 length_l90.py base 1800 results/length_l90.json`, `KF=8 python3 length_l90.py ppm 1200 results/length_l90.json results/length_ppm.json` | `gen_dumps.py`의 세기 이동 S, 길이별 L90 범위, 길이별 합성 시험 페이지 7장 | 약 4분 + 약 6분 |
| `length_h0md.py` | 합성 시험 H0 조건(−120 / −88 / −82 / −40 dBm, 멀티패스)의 MD 규칙 오검출률, 길이별: `KF=8 python3 length_h0md.py 2000 results/length_h0md.json` | 길이별 MD 규칙 통계 범위 | 약 10분 |
| `expect_synthetic.py` | `run_synthetic_test.m`의 파이썬 판: 5열 94개·2열 49개 판정을 같은 이름으로 계산 (세트의 `synthetic_set.txt`에서 길이·세기 이동·이동량을 읽음): `python3 expect_synthetic.py <root5> <root2> results/expect_n64.json` | 길이별 기대값 (`results/expect_n*.json`) | 수십 초 |
| `drift_check.py` | 이동량 0.2496 / 0.6 / 1.2탭에서 이동 폴더(약한 세기 800개, 강한 세기 80개)의 이동 탐색 Pd, 고른 기울기, 기울기, 제안 후보 수, 길이별: `python3 drift_check.py <work dir> results/drift_check.json` | 이동 판정 범위와 확인 구간 | 약 90분 |
| `build_expect.py` | 길이별 판정 범위를 만들어 MATLAB 함수 `test/synthetic_expect.m`으로 씀 (32심볼은 기존 범위 그대로): `python3 build_expect.py results ../test/synthetic_expect.m`, 범위와 파이썬 값 대조는 `python3 build_expect.py results --report` | `run_synthetic_test`의 길이별 범위 | 수 초 |

`results/`에는 위 스크립트의 결과(JSON)를 둡니다. 길이별 기준 세트는 `gen_dumps.py <root> --N <N>`과 `make_2col.py`로 만들고(5열 약 100 MB, 2열 약 60 MB씩, 32심볼 약 7분 ~ 256심볼 약 70분), 저장소에는 넣지 않습니다.

dBm은 MD 문서의 환산(−100 dBm ↔ 펄스당 SNR −11.4 dB, 즉 A²[dB] = dBm + 82.5)을 쓴 모델 값입니다.
