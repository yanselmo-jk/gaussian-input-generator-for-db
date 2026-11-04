# rmlipdata

Li/Na 배터리 전해질용 reactive MLIP(machine-learned interatomic potential) 학습 data를 만드는 도구.

- `rmd generate`: GFN2-xTB(gas phase)로 구조와 반응 경로를 만들고, 각 구조를 Gaussian16 입력(`.gjf`)으로 씀.
  - 계산 수준: wB97X-D/def2-TZVPPD
  - 함께 받는 값: force, Mulliken / NBO(NPA) / Loewdin charge
- `rmd parse`: Gaussian16 log에서 atomic position, energy, force, charge, spin을 읽어 extxyz로 저장.
- 생성한 모든 구조의 id와 metadata는 SQLite 파일 하나(`registry.sqlite`)에 기록됨.

## 설치

```bash
conda create -n rmlipdata -c conda-forge python=3.11 "xtb>=6.7" ase rdkit numpy pyyaml pytest
conda activate rmlipdata
pip install -e .
pytest -q tests        # 7 passed
```

## 사용 흐름

```bash
rmd generate configs/example.yaml --dry   # xtb 없이 sampler별 대상 수만 출력
rmd generate configs/example.yaml         # 구조 생성, .gjf와 registry 작성
# runs/example/manifest.csv의 .gjf를 클러스터에 제출 (제출 스크립트는 각자 환경 것 사용)
rmd parse logs/ -o data.extxyz --registry runs/example/registry.sqlite
rmd status --registry runs/example/registry.sqlite
```

- `generate`는 같은 config로 다시 실행해도 이미 있는 구조를 다시 만들지 않음 (id로 중복 제거).
- `--only bond_scan reaction`처럼 sampler 일부만 실행할 수 있음.
- 이미 만든, 아직 계산하지 않은(`generated`) deck의 `%nprocshared`, `%mem`을 바꾸려면 아래 명령을 씀. `.gjf`와 `manifest.csv`가 함께 바뀜.
  ```bash
  rmd resources --registry runs/example/registry.sqlite --nproc 32 --mem 64GB
  rmd resources --registry runs/example/registry.sqlite --nproc 8 --mem 16GB --sampler monomer pair_scan
  ```

### 산출물 (`out` 폴더, 예: `runs/example/`)

| 파일 | 내용 |
|---|---|
| `registry.sqlite` | 구조 목록과 metadata (아래 "Registry" 참고) |
| `manifest.csv` | 아직 계산하지 않은(`generated`) `.gjf` 목록. 제출 스크립트에 넘기는 용도 |
| `gjf/<sampler>/<id>.gjf` | Gaussian16 입력. 파일명이 곧 registry id |
| `species/` | 분자별 xTB 최적화 구조 cache |

log 파일명은 `.gjf`와 같은 이름(`<id>.log`)을 유지해야 함. parser는 파일명으로 registry의 행을 찾음.

## Config 요약 (`configs/example.yaml`)

| key | 의미 |
|---|---|
| `species` | 사용할 양이온(Li, Na), 음이온(PF6, FSI, ...), 용매(EC, DMC, ...). 목록은 `rmlipdata/data/*.tsv` |
| `charge.delta_q` | 전자 수 변화. 0 = 형식 전하 그대로, -1 = 환원(전자 1개 추가), +1 = 산화 |
| `charge.mult_offsets` | 최저 multiplicity에 더할 값. `[0, 2]` = singlet+triplet 또는 doublet+quartet |
| `cluster.max_components` | cluster 하나에 들어가는 분자 + 음이온의 최대 개수 |
| `cluster.max_cations` | cluster 하나에 들어가는 Li+/Na+의 최대 개수 |
| `cluster.max_combos` | 조합 중 무작위로 고를 개수 |
| `limits.max_atoms` | 원자 수 상한 |
| `g16.nproc`, `g16.mem` | Gaussian 입력의 `%nprocshared`, `%mem` 값과 manifest의 같은 열. 새로 생성하는 deck에만 적용 (이미 만든 deck은 `rmd resources`로 변경) |
| `samplers.<이름>.on` | sampler 켜기/끄기 |
| `samplers.<이름>.max_structures` | 그 sampler가 만들 구조 수 상한 (spin 상태로 늘리기 전 기준) |

Sampler 종류:

| sampler | 만드는 구조 | Gaussian job |
|---|---|---|
| `monomer` | 분자/음이온의 최적화 구조, normal-mode 변형, torsion/angle 변형, 단원자 양이온 | opt / sp |
| `pair_scan` | 두 성분 사이 거리를 0.75–2.0배로 바꾼 binding curve | sp |
| `cluster` | M+(solvent)n(anion)m cluster 최적화 구조와 무작위로 흔든 구조 | opt / sp |
| `bond_scan` | bond 하나를 끊어질 때까지 늘인 xTB relaxed scan | sp (BS singlet 포함) |
| `reaction` | 첨가(add), 치환(substitute), F/H 이동(transfer) 경로의 최고점 ±k개 frame (pseudo-TS) | sp (BS singlet 포함) |
| `hop` | Li+/Na+가 한 donor 원자에서 다른 donor 원자로 옮겨 가는 경로 | sp |
| `nanoreactor` | 고온 xTB MD snapshot (기본 off) | sp |

각 구조는 `delta_q`의 각 값과 `mult_offsets`의 각 값으로 늘어남. BS(broken-symmetry) singlet은 `guess=mix`를 씀.

## Parser 출력 (extxyz)

| 구분 | key | 단위 / 의미 |
|---|---|---|
| 좌표 | `positions` | Å |
| energy | `energy` | eV (SCF energy) |
| force | `forces` | eV/Å |
| charge | `mulliken_charges`, `lowdin_charges`, `nbo_charges` | e. 해당 frame에 출력이 없으면 NaN |
| spin | `mulliken_spins`, `lowdin_spins`, `nbo_spins` | α−β spin density. restricted closed-shell은 0 |
| info | `id`, `step`, `frame` | registry id, Link1 step, frame 번호 |
| info | `charge`, `spin` | 전체 전하, multiplicity (`spin`은 multiplicity를 뜻함) |
| info | `job`, `converged`, `normal_termination` | job 종류와 수렴/종료 상태 |
| info | `s2`, `s2_dev` | ⟨S²⟩, ⟨S²⟩ − S(S+1) (spin contamination 확인용, 걸러내지는 않음) |
| info | `q_sum_ok`, `spin_sum_ok` | Mulliken charge 합과 spin 합 검사 결과 |

opt job은 최적화의 모든 step이 frame이 됨. Gaussian은 charge를 첫 step과 수렴 후에만 출력하므로 중간 step의 charge는 NaN임.

## Registry (SQLite) 사용법

### SQLite란
SQLite는 표(table)를 파일 하나에 담는 database임. 서버가 필요 없고, `registry.sqlite` 파일을 복사하면 그대로 옮겨짐. 이 파일 안에는 `structures`라는 표가 하나 있고, 한 행이 Gaussian 입력 하나에 해당함.

| column | 의미 |
|---|---|
| `id` | 구조 id (= `.gjf`/`.log` 파일명). 같은 원자 종류, 좌표(1e-3 Å 반올림), charge, mult, BS 여부, job이면 같은 id |
| `parent_id` | 부모 구조 id (현재 미사용) |
| `sampler`, `composition`, `formula`, `n_atoms` | 만든 sampler, 구성(예: `Li1 PF6_1 EC2`), 화학식, 원자 수 |
| `charge`, `mult`, `bs`, `job` | 전하, multiplicity, BS singlet 여부(0/1), `sp` 또는 `opt` |
| `params` | sampler별 상세 정보 (JSON 문자열: scan 좌표, 경로 index, template 등) |
| `gjf_path`, `log_path` | 입력/log 파일 경로 |
| `status` | `generated` (입력만 있음) → `parsed` (log 읽음) 또는 `failed` (비정상 종료) |
| `created` | 생성 시각 |

### 보는 방법
1. 가장 쉬운 방법은 CSV로 내보내 Excel 등으로 여는 것.
   ```bash
   rmd export --registry runs/example/registry.sqlite -o registry.csv
   ```
2. 터미널에서 직접 query.
   ```bash
   sqlite3 runs/example/registry.sqlite
   sqlite> .headers on
   sqlite> SELECT sampler, status, COUNT(*) FROM structures GROUP BY sampler, status;
   sqlite> SELECT id, composition, charge, mult FROM structures WHERE sampler='reaction' AND bs=1 LIMIT 5;
   sqlite> .quit
   ```
3. GUI 프로그램 [DB Browser for SQLite](https://sqlitebrowser.org/)로 파일을 열면 표를 마우스로 볼 수 있음.

자주 쓰는 query:
```sql
-- 계산 실패한 구조
SELECT id, log_path FROM structures WHERE status='failed';
-- 구성별, charge별, multiplicity별 개수
SELECT composition, charge, mult, COUNT(*) FROM structures GROUP BY composition, charge, mult;
-- params(JSON) 안의 값 꺼내기: reaction의 pseudo-TS 최고점만
SELECT id, json_extract(params, '$.template') FROM structures
WHERE sampler='reaction' AND json_extract(params, '$.is_peak')=1;
```

### Python에서 extxyz와 합치기 (raw data DB를 만들 때)
extxyz frame의 `info["id"]`가 registry의 `id`와 같으므로 이것으로 연결함.
```python
import json, sqlite3
from ase.io import read

con = sqlite3.connect("runs/example/registry.sqlite")
con.row_factory = sqlite3.Row
meta = {r["id"]: dict(r) for r in con.execute("SELECT * FROM structures")}

for atoms in read("data.extxyz", ":"):
    row = meta[atoms.info["id"]]
    params = json.loads(row["params"])
    print(row["sampler"], row["composition"], atoms.info["frame"], atoms.get_potential_energy(), params)
```

## 미검증 사항 (test deck으로 확인 예정)
`test_decks/`의 Gaussian 입력 3개를 실행해 log를 받으면 확인하고 parser를 확정함.
1. `iop(6/80=1)`이 Loewdin charge/spin을 출력하는지와 그 출력 형식.
2. open-shell NPA의 alpha/beta block 형식. NBO spin은 α − β natural population으로 계산 중.
3. opt 중간 step의 population 출력 여부.

## 주의
- 이번 xTB(GFN2) 결과에서 Li+–EC의 Li–O 거리는 1.66 Å였음. 문헌의 DFT 값(약 1.75–1.80 Å, 우리 계산 아님)보다 짧음. `test_decks/t1` opt 결과로 확인할 예정.
  - xTB 경로에서 뽑은 SP frame은 이 기하 편향을 그대로 가짐.
  - 단, energy/force label은 DFT 값임.
- 중복 제거는 같은 xtb version과 thread 수에서 결과가 재현된다는 가정에 기댐.
- basis: `rmlipdata/data/basis/def2-TZVPPD.gbs`는 Basis Set Exchange(0.12)에서 받은 것. H, Li, B, C, N, O, F, Na, Si, P, S, Cl 포함.
