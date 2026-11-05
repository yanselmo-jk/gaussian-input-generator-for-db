# rmlipdata

Li/Na 배터리 전해질용 reactive MLIP 학습 data를 만드는 도구.
GFN2-xTB로 구조와 반응 경로를 만들고, 각 구조를 Gaussian16 입력(wB97X-D/def2-TZVPPD)으로 씀.
계산이 끝난 log에서 energy, force, charge, spin을 읽어 extxyz로 저장함.

## 설치

```bash
conda create -n rmlipdata -c conda-forge python=3.11 "xtb>=6.7"
conda activate rmlipdata
python -m pip install -e .      # numpy, ase, rdkit, pyyaml, pytest가 함께 설치됨
python -m pytest -q tests
```

xtb는 pip로 설치되지 않으므로 conda로 설치함. `pip`, `pytest` 대신 `python -m pip`, `python -m pytest`를 쓰면 다른 환경(base 등)의 pip/pytest가 실행되는 일을 막을 수 있음.

## 사용법

```bash
rmd generate configs/example.yaml --dry     # 1. 생성될 대상 수 확인 (xtb 실행 안 함)
rmd generate configs/example.yaml           # 2. 구조 생성 → runs/example/gjf/*.gjf
#                                             3. runs/example/manifest.csv의 .gjf를 클러스터에 제출
rmd submit runs/example/manifest.csv --registry runs/example/registry.sqlite   #    제출한 deck 표시 (manifest에서 빠짐)
rmd parse logs/ -o data.extxyz --registry runs/example/registry.sqlite   # 4. log → extxyz
rmd status --registry runs/example/registry.sqlite                       # 5. 진행 상황
```

log 파일명은 `.gjf`와 같게 유지해야 함 (`<id>.gjf` → `<id>.log`).

## 명령어

| 명령 | 하는 일 |
|---|---|
| `rmd generate CONFIG [--dry] [--only SAMPLER ...]` | 구조 생성, Gaussian 입력 작성 |
| `rmd parse PATH... -o OUT.extxyz [--registry DB] [--keep-all]` | Gaussian log를 extxyz로 변환. OMol25 기준 force/⟨S²⟩ 필터 적용 |
| `rmd submit TARGET... --registry DB [--undo]` | 제출한 deck 표시 (generated → submitted). TARGET: id, .gjf, 폴더, manifest .csv |
| `rmd list --registry DB [--status S] [--sampler ...]` | 상태별 deck 목록 (기본 submitted) |
| `rmd status --registry DB` | sampler별, 상태별 개수 |
| `rmd export --registry DB -o FILE.csv` | registry를 CSV로 내보내기 |
| `rmd resources --registry DB [--nproc N] [--mem 64GB] [--sampler ...]` | 이미 만든 입력의 `%nprocshared`/`%mem` 변경 |
| `rmd species [--kind ...] [--group ...]` | 지원 분자 목록 |
| `rmd view NAME [--opt] [--save FILE]` | 분자 구조를 ASE GUI로 보기 |
| `rmd view DECK.gjf... \| DIR [--save FILE]` | 생성된 Gaussian 입력의 구조를 ASE GUI로 보기 (여러 개면 frame으로 넘겨 봄) |
| `rmd add-species FILE --name NAME [--charge Q]` | `.mol`/`.mol2`/`.pdb`/`.xyz`로 분자 추가 |

## 상세 설명

- [docs/generation.md](docs/generation.md): config, sampler, charge/multiplicity, Gaussian 입력
- [docs/species.md](docs/species.md): 지원 분자, 분자 추가, 구조 보기 (분자, 생성된 `.gjf`)
- [docs/parsing.md](docs/parsing.md): log parser와 extxyz 형식
- [docs/registry.md](docs/registry.md): SQLite registry 사용법 (DB 입문)
