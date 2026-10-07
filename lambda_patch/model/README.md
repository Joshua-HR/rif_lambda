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
| `gen_dumps.py` | 합성 5열 덤프 생성: `python3 gen_dumps.py <root>` (`<root> h1_drift`: 그 폴더만 추가) | `synthetic_dumps/` | 약 6분 |
| `verify_dumps.py` | MATLAB 패치 계산의 파이썬 판: `python3 verify_dumps.py <root>` | `run_synthetic_test` 기대값 | 수 초 |
| `make_2col.py` | 5열 세트에서 같은 CIR의 2열 세트 만들기: `python3 make_2col.py ../../synthetic_dumps ../../synthetic_dumps_2col` | `synthetic_dumps_2col/` | 수 초 |
| `web_results.py` | 5열·2열 세트 분석 결과 JSON (`export_web_results.m`과 같은 형식): `python3 web_results.py ../../synthetic_dumps ../../synthetic_dumps_2col out.json` | 웹 10장 "1단계 결과" 기본 데이터 | 약 10초 |
| `packet.py` | `rif_packet.m`과 `check_packet_assumptions.m`의 파이썬 판 (패킷 판정, 이동 기울기, AGC, 독립성) | `run_synthetic_test` 10단계 기대값 | 수 초 |
| `medium_level.py` | FiRa medium(패킷 10⁻⁶) 결합 규칙별 패킷 Pd와 H0 근거: `KF=8 python3 medium_level.py h1 4000 -122 -106 out.json`, `... h0 4000 out.json` | 설계 문서 9.3 첫 표 | 약 3분 + 1분 |
| `policy_k8.py` | 8개 RIF 정책 비교 (strict, soft, 일관성 하한): `python3 policy_k8.py 3000 -121 -104 out.json` | 설계 문서 9.3 정책 표 | 약 2분 |
| `single_len.py` | RIF 1개 설정, 길이별 L90 (1e-3, 10⁻⁶): `python3 single_len.py 3000 out.json` | 설계 문서 9.3 마지막 표 | 약 1분 |

dBm은 MD 문서의 환산(−100 dBm ↔ 펄스당 SNR −11.4 dB, 즉 A²[dB] = dBm + 82.5)을 쓴 모델 값입니다.
