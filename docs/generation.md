# 구조 생성 (`rmd generate`)

## 흐름

1. config의 `species`를 조합해 cluster 후보를 만듦.
2. 켜진 sampler마다 GFN2-xTB(gas phase)로 구조나 경로를 만듦.
3. 구조마다 charge와 multiplicity 조합으로 늘림.
4. 조합마다 registry에 id를 등록하고, 새 id일 때만 `.gjf`를 씀.
5. 마지막에 `manifest.csv`를 다시 씀.

같은 config로 다시 실행하면 이미 있는 구조는 건너뜀. `--only bond_scan reaction`처럼 sampler 일부만 실행해도 cluster 구조는 같게 나옴 (cluster마다 seed가 고정됨).

## 산출물 (`out` 폴더)

| 파일 | 내용 |
|---|---|
| `registry.sqlite` | 구조 목록과 metadata ([registry.md](registry.md)) |
| `manifest.csv` | 아직 계산하지 않은(`generated`) 입력 목록. 열: id, gjf_path, sampler, n_atoms, charge, mult, bs, job, nproc, mem |
| `gjf/<sampler>/<id>.gjf` | Gaussian16 입력. 파일명 = registry id |
| `species/` | 분자별 xTB 최적화 구조 cache |

## Config (`configs/example.yaml`)

| key | 의미 | 예시 값 |
|---|---|---|
| `out` | 출력 폴더 | `runs/example` |
| `seed` | 난수 seed | `0` |
| `species.cations / anions / solvents` | 사용할 분자 이름 ([species.md](species.md)) | `[Li, Na]` |
| `charge.delta_q` | 전자 수 변화. 0 = 형식 전하, -1 = 환원(전자 1개 추가), +1 = 산화 | `[0, -1, 1]` |
| `charge.mult_offsets` | 최저 multiplicity에 더할 값 | `[0, 2]` |
| `charge.range` | 허용 전체 전하 범위. 벗어나면 xTB 계산 전에 건너뜀 | `[-2, 2]` |
| `g16.nproc`, `g16.mem` | 입력의 `%nprocshared`, `%mem` (새로 만드는 입력에만 적용) | `16`, `32GB` |
| `limits.max_atoms` | 원자 수 상한 | `60` |
| `xtb.threads`, `xtb.timeout_s` | xtb thread 수와 job 1개 제한 시간 | `4`, `1800` |
| `cluster.min_components`, `max_components` | cluster 하나의 분자 + 음이온 개수 범위 | `1`, `3` |
| `cluster.max_cations` | cluster 하나의 Li+/Na+ 최대 개수 | `2` |
| `cluster.max_combos` | 전체 조합 중 무작위로 고를 개수 | `50` |
| `samplers.<이름>.on` | sampler 켜기/끄기 | `true` |
| `samplers.<이름>.max_structures` | sampler 하나가 만들 구조 수 상한 (charge/spin으로 늘리기 전 기준) | 없음 |

## Sampler

| sampler | 만드는 구조 | 주요 설정 | job |
|---|---|---|---|
| `monomer` | 분자/음이온 최적화 구조, normal-mode 변형, torsion/angle 변형, 단원자 양이온 | `n_normal_mode`, `n_distort`, `T_K` | opt / sp |
| `pair_scan` | 두 성분의 거리를 바꾼 binding curve | `max_pairs`, `n_orient`, `scales` | sp |
| `cluster` | M+(solvent)n(anion)m 최적화 구조와 흔든 구조 | `n_distort`, `rattle_A` | opt / sp |
| `bond_scan` | bond 하나를 끊어질 때까지 늘인 경로 | `n_points`, `max_bonds_per_cluster`, `bs` | sp |
| `reaction` | 첨가 / 치환 / F·H 이동 경로의 pseudo-TS | `templates`, `n_points`, `k_peak`, `max_per_cluster`, `bs` | sp |
| `hop` | Li+/Na+가 donor 원자 사이를 옮겨 가는 경로 | `n_points`, `max_per_cluster` | sp |
| `nanoreactor` | 고온 xTB MD snapshot (기본 off) | `T_K`, `ncyc`, `every_n_frames` | sp |

### 각 sampler의 방법

- **monomer**
  - normal-mode 변형: xTB Hessian의 mode마다 고전 harmonic 분포 N(0, kT/λ)에서 변위를 뽑음.
    - 50 cm⁻¹ 미만 mode와 병진·회전 mode는 제외함.
    - 원자 하나라도 0.5 Å 넘게 움직이면 다시 뽑음.
  - 비물리적인 bare Li²⁺/Na²⁺는 만들지 않음.
- **pair_scan**
  - 두 성분을 접촉 상태로 놓고 xTB로 최적화함. 이때의 최단 원자 간 거리를 접촉 거리라고 함.
  - 접촉하지 않은 중성 쌍(최단 거리 > 3 Å)은 먼저 2.5 Å까지 붙인 뒤 xTB로 최적화함.
  - 두 번째 성분을 질량중심 축을 따라 옮겨, 최단 원자 간 거리가 `scales` × 접촉 거리가 되게 함.
- **cluster 조립**
  - 양이온이 있으면 각 성분의 donor 원자(O, N, F)가 양이온을 향하도록 배치함.
    - 배위 방향은 사면체 방향 중에서 고름.
    - 배위 거리는 Li 2.0 Å, Na 2.35 Å.
  - 양이온이 없으면 액체 밀도 크기의 구 안에 무작위로 채운 뒤, 성분들을 2.5 Å 접촉까지 붙임. 기상 xTB 최적화는 몇 Å 떨어진 중성 분자를 끌어당기지 못함.
  - 조립한 다음 xTB로 최적화함.
  - 최적화 후 흩어진 cluster(다른 성분과 4 Å 안에 원자가 하나도 없는 성분이 있음)는 경고와 함께 버림. 예: 양이온 2개짜리 +3 cluster의 쿨롱 폭발.
- **bond_scan**
  - 분자 안의 bond를 대칭 종류별로 하나씩 고름.
  - 0.85 × r0부터 1.05 × (공유결합 반지름 합) + 1.5 Å까지 xTB constrained scan을 함.
  - 이웃 frame 사이에서 xTB energy가 5 eV 넘게 뛰면 그 지점 이후 frame을 버림 (xTB 전자 상태가 바뀐 것). reaction scan에도 같은 규칙을 씀.
- **reaction (pseudo-TS)**
  - template:
    - `add`: 친핵 원자(O/N/F) → 친전자 원자(carbonyl C, P, S, B)
    - `substitute`: add와 동시에 이탈기(O/F)를 끊음
    - `transfer`: F/H가 양이온이나 다른 분자의 친전자 원자로 이동
  - concerted xTB scan의 최고점 ±`k_peak`개 frame에 첫 frame과 마지막 frame을 더해 씀.
  - 엄밀한 TS가 아님. 반응 PES를 sampling하기 위한 구조임.
- **hop**
  - 양이온 위치를 출발점에서 도착 donor 원자 앞까지 선형으로 옮김.
  - 각 점에서 양이온–출발 donor, 양이온–도착 donor, 두 donor 사이 거리를 constraint로 묶고 나머지를 xTB로 최적화함. xtb `$fix`는 이 cluster들에서 원자를 정확히 고정하지 못함 (고정한 원자가 최대 2 Å 움직임).
  - 양이온이 어떤 원자에 공유결합 반지름 합의 0.75배보다 가까워진 frame은 버림.

## Charge와 multiplicity

- charge = 구성 성분 형식 전하의 합 + `delta_q`.
- 최저 multiplicity는 전자 수가 짝수면 1, 홀수면 2. 여기에 `mult_offsets`를 더함.
- `bs: true`인 sampler는 전자 수가 짝수일 때 BS(broken-symmetry) singlet(`guess=mix`)을 추가함. homolytic 절단 경로를 위한 것임.
- xTB 경로는 charge마다 따로 계산함 (uhf = 전자 수 % 2). 같은 경로를 여러 multiplicity로 계산함.

## Gaussian 입력

```
%nprocshared=16
%mem=32GB
#p {R|U}wB97XD/gen nosymm int=ultrafine scf=(xqc,maxcycle=300) pop=(npa) iop(6/80=1) force|opt=(maxcycles=200) [guess=mix]

<id>

<charge> <mult>
<좌표>

<원소별 def2-TZVPPD block>
```

- `R`은 closed-shell singlet에 씀. multiplicity > 1이거나 BS singlet이면 `U`.
- `nosymm`을 써서 force가 입력 좌표계로 나옴.
- `pop=(npa)`는 NBO natural charge, `iop(6/80=1)`은 Loewdin population을 출력하게 함.
- basis는 `rmlipdata/data/basis/def2-TZVPPD.gbs`에서 구조에 있는 원소만 붙임.
  - 출처: Basis Set Exchange 0.12.
  - 포함 원소: H–Kr(Z 1–36)의 36개 원소. def2에서 이 원소들은 ECP 없이 모든 전자를 계산함.
  - Rb(Z 37) 이상은 ECP(`genecp`)가 필요해서 지원하지 않음. 이런 원소가 있으면 입력을 쓰지 않고 오류를 냄.
- 불가능한 charge/multiplicity 조합(parity 불일치, multiplicity > 전자 수 + 1)도 오류를 냄.
- 이미 만든 입력의 `%nprocshared`/`%mem` 변경은 `rmd resources`로 함. 대상은 status가 `generated`인 입력이고, `manifest.csv`도 함께 바뀜.

## id와 중복 제거

- id = `<sampler>-<sha1 앞 12자리>`.
- hash 대상: 원소, 좌표(1e-3 Å 반올림), charge, multiplicity, BS 여부, job.
- 같은 xtb version과 thread 수에서 결과가 재현된다는 가정에 기댐.

## 주의

- 이번 xTB 결과에서 Li+–EC의 Li–O 거리는 1.66 Å였음. 문헌의 DFT 값(약 1.75–1.80 Å, 우리 계산 아님)보다 짧음. xTB 경로에서 뽑은 SP frame은 이 기하 편향을 그대로 가짐. energy/force label은 DFT 값임.
- `nanoreactor`는 MD가 도중에 중단되면 재시작 checkpoint가 어긋날 수 있음 (아직 고치지 않음).
