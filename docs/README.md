# RIF Λ 판정기 해설 페이지 (GitHub Pages)

`index.html` 한 파일로 된 정적 페이지입니다. 빌드 과정이 없고, 위젯의 시뮬레이션은 브라우저에서 돕니다.
외부로 나가는 요청은 Google Fonts 하나뿐이며, 사내망에서 막혀도 렌더링을 막지 않고 시스템 글꼴(맑은 고딕 등)로 표시됩니다.

## 사내 GitHub(GitHub Enterprise) Pages로 게시

1. 이 `docs/` 폴더(`index.html`, `.nojekyll`)를 저장소의 게시할 브랜치에 커밋합니다.
2. 저장소 **Settings → Pages → Build and deployment**
   - Source: **Deploy from a branch**
   - Branch: 게시할 브랜치(예: `main`), Folder: **`/docs`** → Save
3. 몇 분 뒤 같은 화면에 Pages 주소가 표시됩니다. 주소 형식은 사내 설정에 따라
   `https://pages.<사내 GitHub 호스트>/<owner>/<repo>/` 또는 `https://<사내 GitHub 호스트>/pages/<owner>/<repo>/` 입니다.

참고

- Pages 메뉴가 보이지 않으면 사이트 관리자가 GitHub Enterprise Server에서 Pages를 켜야 합니다.
- 저장소의 `/docs`를 이미 다른 용도로 쓰고 있다면, 이 두 파일만 담은 별도 브랜치(예: `gh-pages`, 파일을 브랜치 최상위에 두고 Folder: `/ (root)`)나 별도 저장소를 쓰세요.
- `.nojekyll`은 Jekyll 처리를 끄기 위한 빈 파일입니다. 지우지 마세요.
- 로컬에서 보려면 `index.html`을 브라우저로 바로 열면 됩니다.

## 내용

1. 무엇을 판정하나 (CIR 위젯) 2. 기존 판정기 (결정 평면 위젯) 3. 왜 깨지나 (꼬리 확률) 4. 핵심 관찰
5. 백색화 Λ (상관값 구름 위젯) 6. 문턱 계산기 7. 검증 결과 8. 고정 칩 대안 (CIR 추정 Λ̂) 9. 패킷 판정 (FiRa medium, RIF 8개 위젯) 10. 적용 순서 (합성 5열·2열 시험 포함)

10장 "1단계 결과: 5열 덤프 대 2열 덤프" 그림은 기본으로 Python 미러 값을 보여 줍니다. MATLAB의 `export_web_results`가 만든 JSON을 `docs/results/matlab_results.json`으로 커밋하면 게시된 페이지가 그 결과를 기본으로 보여 주고, 그림의 버튼으로 Python 미러와 바꿔 볼 수 있습니다.

설계 문서: `../lambda_patch/LAMBDA_PATCH_DESIGN.md`. 페이지의 수치는 `../lambda_patch/model/`의 스크립트로 다시 만들 수 있습니다.
