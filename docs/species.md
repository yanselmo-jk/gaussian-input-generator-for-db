# 분자 목록, 추가, 보기

## 지원 분자

```bash
rmd species                      # 전체
rmd species --kind anion         # cation | anion | solvent
rmd species --group gas          # group으로 거르기
```

| 출처 파일 | 내용 | group |
|---|---|---|
| `rmlipdata/data/anions.tsv` | 전해질 음이온 19종 (PF6, FSI, TFSI, BF4, ClO4, BOB, DFOB, ...). `1_electrolyte_inverse_design` commit `318e8ac`에서 복사. 원본에서 제외 표시된 AsF6⁻, Br⁻은 basis가 있어도 계속 제외함 | `anion` |
| `rmlipdata/data/solvents.tsv` | 용매·첨가제 43종 (EC, PC, DMC, EMC, DEC, FEC, DME, VC, ...). 같은 출처 | 원본의 `family` 열 |
| `rmlipdata/data/extra.tsv` | 이 repo에서 추가한 34종 | 아래 표 |
| 코드 내장 | Li+, Na+ | `metal_cation` |

`extra.tsv`의 group (분자 목록은 일반적인 SEI·전해질 문헌 지식으로 고른 것이고, 우리 계산이 아님):

| group | 분자 |
|---|---|
| `impurity` | H2O, HF |
| `hydrolysis` | PF5, POF3, HPO2F2 |
| `gas` | CO2, CO, C2H4, H2, CH4 |
| `sei_inorganic` | LiF, LiOH, Li2O, Li2CO3, NaF, NaOH, Na2O, Na2CO3 |
| `sei_organic` | LEDC, LMC, LiOMe, LiOEt, NEDC, NMC, NaOMe, NaOEt |
| `radical` | EC_ro (고리가 열린 EC radical anion), MeOCO2, CH3, CH3O |
| `glyme` | G3, G4 |
| `lhce_diluent` | BTFE, TTE |

kind에 따른 처리:

| kind | 처리 방식 |
|---|---|
| `cation` | Li+/Na+ 금속 이온. cluster 중심에 놓이고 `hop` 경로를 만듦 |
| `anion` | 전하를 띤 성분 |
| `solvent` | 중성 성분. LiF처럼 Li가 들어 있는 중성 분자도 여기에 속하며 `max_cations` 개수에 들어가지 않음 |

## 분자 추가 (`rmd add-species`)

```bash
rmd add-species my.mol  --name MyMol
rmd add-species ion.xyz --name MyAnion --charge -1
rmd add-species my.mol  --name MyMol --force      # 같은 이름 덮어쓰기
```

| 형식 | 결합과 charge |
|---|---|
| `.mol`, `.sdf`, `.mol2` | 파일에 적힌 결합과 formal charge를 씀. `--charge`가 파일 값과 다르면 오류 |
| `.xyz`, `.pdb` | RDKit이 결합을 추정함. charge 기본값은 0이고 `--charge`로 지정함. 추정에 실패하면 `.mol`/`.mol2`를 쓰라는 메시지가 나옴 |

검사와 동작:
- 이름은 영문자, 숫자, `_`, `+`, `-`만 쓸 수 있음. 기본 분자와 같은 이름은 쓸 수 없음.
- 모든 원소가 H–Kr(Z 1–36) 범위 안에 있어야 함. 그보다 무거운 원소는 ECP가 필요해서 지원하지 않음.
- 3D 좌표가 있어야 함. H가 빠져 있으면 좌표와 함께 추가하고 경고를 냄.
- 조각이 여러 개면 경고만 내고 허용함.
- kind 기본값은 charge < 0이면 `anion`, 아니면 `solvent`. `--kind`로 바꿀 수 있음.
- Li/Na가 들어간 분자는 결합 추정이 잘 안 되므로 `.mol`/`.mol2`를 권장함.

저장 위치: `$RMLIPDATA_HOME` (기본 `~/.local/share/rmlipdata/`)

| 파일 | 내용 |
|---|---|
| `species/<name>.mol` | 저장된 구조 (H 포함, 3D) |
| `user_species.tsv` | 목록. 열: name, smiles, charge, kind, source, added |

추가한 분자는 group `user`로 나옴. config의 `species`에 이름을 넣으면 다른 분자와 똑같이 쓰임. 구조 생성은 파일의 좌표에서 출발해 xTB로 최적화함. `--force`로 바꾸면 cache 이름이 바뀌므로 예전 구조를 다시 쓰지 않음.

## 구조 보기 (`rmd view`)

```bash
rmd view EC                 # ASE GUI 창 (Linux, X display 필요)
rmd view EC --opt           # xTB 최적화 후 표시, energy 출력
rmd view EC --save ec.xyz   # 창 없이 파일로 저장
```

- 기본 분자는 SMILES에서 만든 초기 구조를 보여 주고, 추가한 분자는 저장된 파일 구조를 보여 줌.
- display가 없으면 안내 메시지가 나옴:
  - WSL2(WSLg): `export DISPLAY=:0`
  - ssh: `ssh -X`
  - 창 없이: `--save`
