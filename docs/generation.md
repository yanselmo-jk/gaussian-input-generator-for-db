# 구조 생성 (`rmd generate`)

## 흐름

1. config의 `species`를 조합해 cluster 후보를 만듦.
2. 켜진 sampler마다 GFN2-xTB(gas phase)로 구조나 경로를 만듦.
   - delta_q와 sampler 한 조합을 stage라 함. stage는 `xtb.n_max_process`개까지 별도 process로 동시에 실행됨.
   - 난수는 seed, sampler 이름, 조합(또는 분자, 쌍), charge로만 정해짐. 동시 실행 수, `--only`, 앞 조합의 xTB 실패와 상관없이 같은 구조가 나옴.
3. 구조마다 charge와 multiplicity 조합으로 늘림.
4. 조합마다 registry에 id가 없거나 deck 경로가 비어 있으면 `.gjf`를 임시 파일로 쓰고 이름을 바꾼 뒤 registry에 등록함. 중간에 멈춰도 deck 없는 행이 남지 않음.
   - 그 전에, 같은 sampler·delta_q stage 안에서 같은 조성·charge·원소 순서의 이전 구조와 Kabsch 정렬 후 원자 최대 변위가 `limits.dedup_A`(기본 0.03 Å) 미만인 구조는 버림. 끝의 요약에 sampler별로 버린 개수가 나옴. 다른 stage와는 비교하지 않음 (실행 순서와 무관한 결과를 위해).
   - vertical deck과 stage 구조가 같은 deck(같은 id)이면 stage 구조의 params가 남고 `parent_id`는 비움.
5. 마지막에 `manifest.csv`를 다시 씀.

example config는 불순물·생성물(H2O, HF, PF5, POF3, CO2, LiF, LEDC, VC, CH3O, EC_ro), Δq −2, sampler별 고스핀 비율(평형 쪽 0.1, scan 1.0, hop 0.25)을 포함함.

실행 중 경고는 `<out>/warnings.log`에 기록되고, 끝에 종류별 개수만 출력됨. stage마다 같은 cluster를 다시 만들어 같은 경고가 반복되므로, 여러 stage에서 나온 같은 문장은 한 stage에서 가장 많이 나온 횟수만 남김. 그래서 서로 다른 stage에서 따로 생긴 같은 문장의 경고는 덜 세어질 수 있음. 같은 config로 다시 실행하면 이미 있는 구조는 건너뜀. `--only bond_scan1 reaction`처럼 sampler 일부만 실행해도 cluster 구조는 같게 나옴 (cluster마다 seed가 고정됨).

xTB 결과는 `<out>/xtb_cache`에 저장됨. key는 xtb 입력 전체(좌표, constraint/scan 입력, charge, uhf, etemp, xtb binary 경로와 수정 시각)의 sha1임. 다시 실행하면 xtb를 돌리지 않고 같은 구조를 받아 registry에서 기존 id로 걸러짐. 결과는 JSON 파일(`<key>.json`)로 저장됨. v0.1.9 이전의 pickle cache는 쓰지 않으므로 처음 실행 때 xTB를 다시 계산함. 실패한 계산은 저장하지 않음. xtb를 다시 돌리려면 이 폴더를 지움.

## 산출물 (`out` 폴더)

| 파일 | 내용 |
|---|---|
| `registry.sqlite` | 구조 목록과 metadata ([registry.md](registry.md)) |
| `manifest.csv` | 아직 계산하지 않은(`generated`) 입력 목록. 열: id, gjf_path, sampler, n_atoms, charge, mult, bs, job, nproc, mem |
| `gjf/<job>/<sampler>/<id>.gjf` | Gaussian16 입력. `<job>`은 `sp`(single point + force) 또는 `opt`(최적화). 파일명 = registry id |
| `species/` | 분자별 xTB 최적화 구조 cache |
| `xtb_cache/` | xTB 결과 cache (key = xtb 입력의 sha1) |

## Config (`configs/example.yaml`)

| key | 의미 | 예시 값 |
|---|---|---|
| `out` | 출력 폴더 | `runs/example` |
| `seed` | 난수 seed | `0` |
| `species.cations / anions / solvents` | 사용할 분자 이름 ([species.md](species.md)) | `[Li, Na]` |
| `charge.delta_q` | 전자 수 변화. 0 = 형식 전하, -1 = 환원(전자 1개 추가), +1 = 산화 | `[0, -1, 1]` |
| `charge.mult_offsets` | 최저 multiplicity에 더할 값 | `[0, 2]` |
| `charge.range` | 허용 전체 전하 범위. 벗어나면 xTB 계산 전에 건너뜀 | `[-2, 2]` |
| `charge.redox_fraction` | Δq ≠ 0 구조 중 deck으로 쓸 비율 (0–1). 구조마다 id hash로 정하므로 다시 실행해도 같은 구조가 뽑힘. xTB 계산 수는 줄지 않음 | `1.0` |
| `charge.high_spin_fraction` | 최저값보다 높은 multiplicity deck을 쓸 비율 (0–1). BS singlet은 해당 없음 | `1.0` |
| `charge.vertical_fraction` | sp 구조 중 같은 좌표로 charge q−1, q+1 deck도 쓸 비율 (id hash로 선택). 형식 전하 ±1 안에서만. q−1은 Li/Na 원자가 있고 결과 charge가 −1 이상일 때만, q+1은 Li/Na만으로 된 구조 제외. 단원자 양이온과 양이온 쌍(좌표가 charge와 무관)은 제외. params `vertical_dq`, `parent_id` | `0.2` |
| `g16.nproc`, `g16.mem` | 입력의 `%nprocshared`, `%mem` (새로 만드는 입력에만 적용) | `16`, `32GB` |
| `g16.pop_always` | opt 입력에 `pop=(npa,always)`를 넣어 매 step population을 출력함. false면 `pop=(npa)` | `true` |
| `limits.max_atoms` | 원자 수 상한 | `60` |
| `limits.dedup_A` | 중복 구조 판정 거리 (원자 최대 변위). 0이면 끔 | `0.03` |
| `samplers.<이름>.high_spin_fraction` | 그 sampler에만 적용하는 `charge.high_spin_fraction` | scan sampler `0.25` |
| `samplers.<이름>.max_atoms` | 그 sampler에만 적용하는 원자 수 상한 (`cluster_md`) | `100` |
| `samplers.<이름>.path_rattle_A` | 경로 frame 두 개 중 하나를 σ로 흔든 구조로 바꿈 (bond_scan1/2, reaction/2, hop). params `rattle_A` | `0.07` |
| `samplers.reaction.peak_rattle_n` | 최고점 frame을 흔든 복사본 수. params `peak_rattle` | `2` |
| `xtb.n_max_process` | 동시에 실행할 stage 수. CPU 사용량은 약 n_max_process × threads | `2` |
| `xtb.threads`, `xtb.timeout_s` | xtb job 1개의 thread 수와 제한 시간 | `2`, `1800` |
| `cluster.min_components`, `max_components` | cluster 하나의 분자 + 음이온 개수 범위. `max_components`는 int 또는 양이온 수별 dict | `1`, `{0: 3, 1: 4, 2: 3}` |
| `cluster.cation_weights` | `max_combos`를 양이온 0/1/2개 조합에 나누는 비율. 양이온 2개 이상 조합은 음이온을 반드시 포함 | `[0.2, 0.7, 0.1]` |
| `cluster.max_cations` | cluster 하나의 Li+/Na+ 최대 개수 | `2` |
| `cluster.max_combos` | 전체 조합 중 무작위로 고를 개수 | `50` |
| `samplers.<이름>.on` | sampler 켜기/끄기 | `true` |
| `samplers.<이름>.max_structures` | sampler 하나가 만들 구조 수 상한 (charge/spin으로 늘리기 전 기준) | 없음 |

## Sampler

| sampler | 만드는 구조 | 주요 설정 | job |
|---|---|---|---|
| `monomer` | 분자/음이온 최적화 구조, normal-mode 변형, torsion/angle 변형, 단원자 양이온 | `n_normal_mode`, `n_distort`, `T_K` | opt / sp |
| `atoms` | 단원자 기준 deck: 설정한 종의 모든 원소, charge −1/0/+1 (H⁺ 제외). delta_q 0 stage에서만, vertical deck 없음 | `high_spin_fraction` (1.0 권장: N 4중항 같은 Hund 바닥 상태) | sp |
| `pair_scan` | 두 성분의 거리를 바꾼 binding curve. 양이온 쌍(Li–Li, Li–Na, Na–Na) 포함 | `max_pairs`, `n_orient`, `scales` | sp |
| `cluster` | M+(solvent)n(anion)m 최적화 구조, 흔든 구조, 늘이거나 줄인 구조 | `n_distort`, `rattle_A` (숫자나 목록, 기본 0.05 Å, xTB 기준 원자당 약 0.1 eV), `dilate`, `m_stretch` | opt / sp |
| `cluster_md` | 양이온 1개 + 용매 n개(+ 음이온)의 xTB MD frame을 RMSD로 골라 씀 | `n_solvent`, `with_anion`, `T_K`, `time_ps`, `equil_ps`, `dump_fs`, `n_frames`, `max_clusters`, `max_atoms` | sp |
| `bond_scan1` | bond 하나를 끊어질 때까지 늘이며 xTB로 relax한 scan | `n_points`, `max_bonds_per_cluster`, `bs` | sp |
| `bond_scan2` | 같은 bond를 IDPP interpolation으로 늘인 경로. 끝점만 xTB 최적화. 고리가 아닌 bond는 더 먼 해리 점 | `n_points`, `max_bonds_per_cluster`, `bs`, `far_A` | sp |
| `reaction` | 첨가 / 치환 / F·H 이동 / radical 결합 / H 빼앗기 / 분자 내 H·F 이동 / 양이온 보조 고리 열림의 pseudo-TS (xTB relaxed scan) | `templates`, `n_points`, `k_peak`, `max_per_cluster`, `bs` | sp |
| `reaction2` | reaction과 같은 template. xTB로 최적화한 생성물 끝점까지 IDPP interpolation | `templates`, `n_points`, `max_per_cluster`, `bs`, `free_opt` | sp |
| `hop` | Li+/Na+가 donor 원자 사이를 옮겨 가는 경로, 첫 배위권 donor가 양이온에서 멀어지는 경로(ligand exchange) | `n_points`, `max_per_cluster`, `exchange`, `exchange_max_per_cluster` | sp |
| `nanoreactor` | 고온 xTB MD snapshot (기본 off) | `T_K`, `ncyc`, `every_n_frames` | sp |

### sampler별로 만드는 deck

xTB 구조 하나는 charge와 spin 조합 수만큼 deck(Gaussian 입력 파일)으로 늘어남.
- **charge**: `charge.delta_q` 값마다 "형식 전하 + Δq"로 하나씩 만듦. 결과가 `charge.range` 밖이면 xTB 계산 전에 건너뜀.
- **multiplicity**: 최저 multiplicity에 `charge.mult_offsets`를 더한 값들을 씀. 기본은 최저값과 최저값+2임. 전자 수 + 1보다 큰 multiplicity는 뺌.
- **BS singlet**: `bs: true`인 sampler(`bond_scan1`, `bond_scan2`, `reaction`, `reaction2`)에서 전자 수가 짝수인 구조의 singlet은 R 대신 U + `guess=mix`(BS singlet)로만 만듦. 같은 구조에 singlet label이 하나만 생기게 하기 위함임. OMol25도 bond가 끊어질 수 있는 singlet을 UKS와 대칭을 깬 초기 guess로 계산함. bond가 끊어지지 않은 구조에서는 UKS 해가 R 해로 돌아옴.
- **상한**: `max_structures`는 charge와 spin으로 늘리기 전의 구조 수에 대한 상한임.

| sampler | opt deck | sp deck | BS singlet |
|---|---|---|---|
| `monomer` | 분자와 음이온마다 xTB 최적화 구조 1개. Gaussian이 다시 최적화하며, 최적화의 모든 step이 frame이 됨 | normal-mode 변형 `n_normal_mode`개와 torsion/angle 변형 `n_distort`개. Li+/Na+는 원자 1개짜리 sp만 만들고, 2+ 상태는 만들지 않음 | 없음 |
| `pair_scan` | 없음 | 쌍 × 방향(`n_orient`) × 거리(`scales`)마다 1개 | 없음 |
| `cluster` | 조합마다 조립 후 xTB로 최적화한 cluster 1개 | 흔든 구조 `rattle_A` 값마다 `n_distort`개, 늘이거나 줄인 구조 `dilate` 값마다 1개 | 없음 |
| `cluster_md` | 없음 | 조성 × 온도마다 MD frame에서 최대 `n_frames`개 | 없음 |
| `bond_scan1` | 없음 | bond마다 scan frame 최대 `n_points`개. energy가 급변하거나 cluster가 분리된 frame부터 뒤는 버림 | `bs: true`면 singlet을 BS로만 |
| `bond_scan2` | 없음 | bond마다 최대 `n_points` − 1개. 첫 image는 cluster와 같아서 빼고, 원자가 겹치거나 다른 bond가 1.3배 넘게 늘어난 image도 버림. 고리가 아닌 bond는 `far_A` 거리마다 1개 추가 | `bs: true`면 singlet을 BS로만 |
| `reaction` | 없음 | 경로마다 끝 frame과 최고점 ±`k_peak` frame (시작 frame은 이완된 cluster와 같아 뺌) | `bs: true`면 singlet을 BS로만 |
| `reaction2` | 없음 | 경로마다 최대 `n_points` − 1개. 첫 image는 빼고, 원자가 겹치거나 절단 대상이 아닌 bond가 1.3배 넘게 늘어난 image는 버림 | `bs: true`면 singlet을 BS로만 |
| `hop` | 없음 | 경로마다 최대 `n_points` − 1개 (t = 0은 이완된 cluster라 뺌). 양이온이 다른 원자와 겹치는 점은 버림. exchange 경로마다 첫 점을 뺀 scan frame | 없음 |
| `nanoreactor` | 없음 | MD trajectory에서 `every_n_frames`마다 snapshot 1개 | 없음 |

opt deck은 `g16.pop_always: true`(기본값)일 때 매 step의 charge와 spin을 출력함. 자세한 내용은 [parsing.md](parsing.md)에 있음.

### 각 sampler의 방법

- **monomer**
  - normal-mode 변형: xTB Hessian의 mode마다 고전 harmonic 분포 N(0, kT/λ)에서 변위를 뽑음.
    - 50 cm⁻¹ 미만 mode와 병진·회전 mode는 제외함.
    - 원자 하나라도 0.5 Å 넘게 움직이면 다시 뽑음.
  - 비물리적인 bare Li²⁺/Na²⁺는 만들지 않음.
  - 형식 전하가 아닌 q에서 xTB 최적화 후 공유결합 조각 수가 형식 전하 구조와 다르거나 조각(Li/Na 포함)이 4 Å 넘게 떨어지면 경고하고 건너뜀.
  - normal-mode와 torsion/angle 변형은 분자·charge마다 따로 정한 난수를 씀 (Hessian 실패가 변형 구조를 바꾸지 않음).
- **pair_scan**
  - 두 성분을 접촉 상태로 놓고 xTB로 최적화함. 이때의 최단 원자 간 거리를 접촉 거리라고 함.
  - 접촉하지 않은 중성 쌍(최단 거리 > 3 Å)은 먼저 2.5 Å까지 붙인 뒤 xTB로 최적화함.
  - 두 번째 성분을 질량중심 축을 따라 옮겨, 최단 원자 간 거리가 `scales` × 접촉 거리가 되게 함.
  - 최적화 후 두 성분 사이에 공유결합(거리 < 1.2 × 공유결합 반지름 합, Li/Na 쌍 제외)이 생긴 방향은 scan하지 않고 경고를 남김. 결합한 adduct를 통째로 옮기면 새 bond 길이에 scale이 곱해지기 때문임 (예: EC + CO₂ q=−1에서 접촉 1.31 Å, scale 0.75면 0.98 Å).
  - 양이온 쌍은 xTB 없이 공유결합 반지름 합(가정, Li–Li 2.56 Å)을 접촉 거리로 쓰고 방향은 1개만 만듦. params `contact: rcov_sum`.
  - 양이온만으로 된 계(양이온 쌍, 양이온만 있는 cluster)는 Δq > 0(추가 산화)을 만들지 않음.
- **cluster 조립**
  - 양이온이 있으면 각 성분의 donor 원자(O, N, F)가 양이온을 향하도록 배치함.
    - 배위 방향은 사면체 방향 중에서 고름.
    - 배위 거리는 Li 2.0 Å, Na 2.35 Å.
  - 양이온이 없으면 액체 밀도 크기의 구 안에 무작위로 채운 뒤, 성분들을 2.5 Å 접촉까지 붙임. 기상 xTB 최적화는 몇 Å 떨어진 중성 분자를 끌어당기지 못함.
  - 조립한 다음 xTB로 최적화함.
  - xTB SCC가 수렴하지 않으면(기체상 radical 음이온, dianion 등) `--etemp 1000`(Fermi smearing)으로 한 번 다시 돌림. 그 결과를 쓴 record의 params에 `xtb_etemp_K: 1000`이 남음. smearing은 xTB 구조 생성에만 쓰이고 DFT label과는 무관함. 이 재시도는 모든 xTB 계산에 적용됨.
  - 다시 돌린 opt에서 성분이 4 Å 넘게 떨어지면 opt 전 접촉 구조를 씀 (기체상에서 결합 구조가 없는 계, 예: EC_ro⁻ + CO₂•⁻).
  - 최적화 후 흩어진 cluster(다른 성분과 4 Å 안에 원자가 하나도 없는 성분이 있음)는 경고와 함께 버림. 예: 양이온 2개짜리 +3 cluster의 쿨롱 폭발.
- **cluster**
  - `rattle_A` 목록: σ마다 `n_distort`개. params `sigma_A`.
  - `dilate`: 성분을 강체로 두고 질량중심만 cluster 중심에서 배율만큼 옮김. params `kind: dilate`, `factor`. 압축(0.9)은 반발 영역, 팽창(1.1, 1.25)은 해리 쪽 구조.
  - `m_stretch`: 양이온 1개인 cluster에서 리간드를 강체로 양이온에서 멀어지게 옮겨 최단 양이온–리간드 거리를 배율만큼 늘림. GFN2-xTB의 M–O 거리가 짧기 때문 (Li–O 1.66 Å, Na–O 1.89 Å, xTB, 우리 계산). params `kind: m_stretch`, `factor`.
- **cluster_md**
  - 조성: 양이온 1개 + 용매 `n_solvent`개. `n_solvent`는 목록 또는 양이온별 dict (`{Li: [3, 4, 5], Na: [4, 5, 6]}`). `with_anion: true`면 용매 하나를 음이온으로 바꾼 조성도 만듦. 원자 수는 `max_atoms`(없으면 `limits.max_atoms`) 이하만.
  - 조립(cluster와 같은 배치) 후 원점 중심으로 옮겨 xTB MD(GFN2, `T_K`마다 `time_ps`). 구 벽(wall) 반지름 = 원점에서 가장 먼 원자까지 거리 + 2 Å, 성분 증발 방지.
  - `equil_ps` 이전 frame은 버림. 성분이나 공유결합 조각(frame마다 다시 구함)이 다른 원자에서 4 Å 넘게 떨어진 frame, 원자가 0.6 × (공유결합 반지름 합)보다 가까운 frame도 버림.
  - 남은 frame에서 RMSD max-min(이미 고른 frame과 RMSD가 가장 먼 frame을 차례로)으로 `n_frames`개. OMol25 전해질 샘플링(MD + RMSD 다양성 선택)과 같은 방식이지만 OPLS4 대신 xTB MD임.
  - MD가 실패하면 경고를 남기고 건너뜀. params `T_K`, `frame`, `rmsd_min_A`, `names`.
  - 비용: 조성 하나에 1 ps당 약 27 s (c4 test config, xtb thread 1개, 2026-11-05, 우리 계산). MD가 generate 시간 대부분을 차지함.
- **bond_scan1**
  - 분자 안의 bond를 대칭 종류별로 하나씩 고름.
  - 0.85 × r0부터 1.05 × (공유결합 반지름 합) + 1.5 Å까지 xTB constrained scan을 함.
  - 이웃 frame 사이에서 xTB energy가 5 eV 넘게 뛰면 그 지점 이후 frame을 버림. reaction scan에도 같은 규칙을 씀.
    - example config 기준 scan 387개 중 2개(0.5%)에서 일어남 (2026-11-05, 우리 계산).
    - 원인 1: scan 최적화가 멈춰 같은 좌표에 잘못된 energy가 기록됨. 같은 좌표를 새로 계산하면 jump가 없음.
    - 원인 2: 한 step 안에서 다른 energy 바닥으로 넘어감. 예: ClO4⁻의 Cl–O 절단 중 떨어진 O가 용매의 H를 빼앗음 (−6.45 eV, 새로 계산해도 같음).
  - relaxed scan이라 대상 bond 외의 원자는 각 점에서 xTB로 이완됨. 성분 간 거리가 4 Å를 넘는 첫 frame부터 버림 (예: q=+2 Li⁺–EC•⁺에서 C=O를 늘이면 Li⁺가 밀려남). reaction scan에도 같은 규칙을 씀.
- **bond_scan2**
  - bond는 bond_scan1과 같은 규칙으로 고름.
  - 끝점: bond의 j쪽 조각(고리 bond면 j 원자 하나)을 bond 방향으로 옮겨 i–j 거리를 r1 = 1.05 × (공유결합 반지름 합) + 1.5 Å로 만든 뒤, 그 거리를 고정하고 xTB로 최적화함.
  - cluster와 끝점 사이를 IDPP interpolation(`ase.mep`)으로 `n_points`개 image로 이음. 중간 image는 energy model 없이 기하로만 정해지므로 xTB scan의 jump나 재배열이 없음. OMol25의 geodesic / Popcornn 경로와 같은 방식임.
  - 첫 image(cluster 자체)는 cluster opt deck과 같아서 쓰지 않음. 두 원자가 0.7 × (공유결합 반지름 합)보다 가까운 image는 버림. 끝점 최적화에서 성분이 4 Å 넘게 떨어지면 그 bond를 건너뜀.
  - 끝점에서 대상 bond 외의 bond가 바뀌면(재배열) 고리가 아닌 bond는 최적화 전 구조(조각만 옮긴 구조)를 끝점으로 씀(`endpoint: rigid`). 그 구조에서도 새 bond가 생기면 건너뜀. 고리 bond는 건너뜀. c4 test config에서 xTB 끝점 16개 중 8개가 H 빼앗기, F 추가 이탈 같은 실제 반응을 일으킴 (2026-11-05, 우리 계산).
  - 다른 bond가 원래 길이의 1.3배를 넘는 image도 버림.
  - `far_A`: 고리가 아닌 bond는 끝점의 j쪽 조각을 더 옮겨 i–j가 `far_A`인 구조도 씀 (완전 해리 쪽). 원자가 겹치거나, 다른 bond가 1.3배 넘게 늘어나거나, 새 bond가 생기면 버림. params `method: far`, `far: true`.
  - params: `method: idpp`, `r_A`, `index`, `ring`, `endpoint`, `E_xtb_end_eV`(xTB 끝점일 때만).
- **reaction (pseudo-TS)**
  - template:
    - `add`: 친핵 원자(O/N/F) → 친전자 원자(carbonyl C, P, S, B)
    - `substitute`: add와 동시에 이탈기(O/F)를 끊음
    - `transfer`: F/H가 양이온이나 다른 분자의 친전자 원자로 이동
    - `radical_couple`: 서로 다른 성분의 radical 원자끼리 결합
    - `h_abstract`: radical 원자가 다른 성분의 H를 빼앗음
    - `intra_shift`: 한 분자 안에서 H/F가 2–3 bond 떨어진 원자로 이동
    - `ring_open_M`: 양이온이 고리 O(이중결합 없는 O)에 다가가며 고리 C–O가 끊어짐
    - `h_transfer_ox`: 산화(Δq > 0) cluster에서만. 한 성분의 C–H의 H가 다른 중성 성분의 carbonyl O나 음이온의 F/O/N으로 이동. 출처: Borodin & Xing, PCCP 14, 12838 (2012) (EC 라디칼 양이온의 ethylene H를 이웃 EC carbonyl O가 가져감); Xing et al., J. Phys. Chem. A 115, 13896 (2011) (PC–PF6 산화 첫 단계에서 HF + PF5)
  - radical 원자는 SMILES의 radical 표시로 찾음 (현재 EC_ro, MeOCO2, CH3, CH3O).
  - Δq ≠ 0이면 xTB 부분전하 변화로도 찾음: q와 q − Δq에서 xTB sp를 돌려, 환원이면 전하가 0.15 e 이상 줄어든 원자, 산화면 0.15 e 이상 늘어난 원자 (H, Li, Na 제외, 최대 2개). SMILES에 radical 표시가 없는 radical 이온(예: 환원된 EC)도 `radical_couple`, `h_abstract`에 쓰임. params `charge_radicals`.
  - 후보는 template별로 섞고 template 순서도 무작위로 정한 뒤 돌아가며 `max_per_cluster`개를 고름. `intra_shift` 후보가 다른 template보다 훨씬 많아 다른 template이 뽑히지 않던 문제 때문임.
  - concerted xTB scan의 최고점 ±`k_peak`개 frame에 첫 frame과 마지막 frame을 더해 씀.
  - 엄밀한 TS가 아님. 반응 PES를 sampling하기 위한 구조임.
- **reaction2**
  - 후보와 선택 규칙은 reaction과 같음.
  - 끝점: 모든 constraint 거리를 목표값으로 고정하고 xTB로 최적화함. 다음 중 하나면 경고와 함께 건너뜀.
    - 출발 cluster가 이미 생성물 상태임.
    - constraint 거리가 목표에서 0.3 Å 넘게 벗어남.
    - 만들어야 할 bond가 생기지 않았거나 끊어야 할 bond가 남음.
    - 그 밖의 bond가 바뀜(재배열).
  - `free_opt: true`면 제약 없이 한 번 더 최적화하고, bond 구성이 같으면 그 구조를 끝점으로 씀(`endpoint: free`). 아니면 `constrained`.
  - cluster와 끝점 사이를 IDPP로 `n_points`개 image로 잇고, energy 없이 기하로 만든 경로라 peak를 고르지 않고 모든 image를 씀.
  - params: `template`, `atoms`, `pairs`, `index`, `endpoint`, `E_xtb_end_eV`.
- **hop**
  - 양이온 위치를 출발점에서 도착 donor 원자 앞까지 선형으로 옮김.
  - 도착 donor는 출발 donor와 다른 원자만 고름 (양이온이 가장 가까운 donor에서 2.6 Å보다 멀면 그 donor가 후보 범위에 들어오기 때문).
  - 각 점에서 양이온–출발 donor, 양이온–도착 donor, 두 donor 사이 거리를 constraint로 묶고 나머지를 xTB로 최적화함. xtb `$fix`는 이 cluster들에서 원자를 정확히 고정하지 못함 (고정한 원자가 최대 2 Å 움직임).
  - 기준 원자 r3(경로 중점 6 Å 이내 heavy atom 중 donor 축에서 가장 먼 원자, 축에서 0.5 Å 이상)와 양이온·두 donor 사이 거리도 고정함 (총 6개). 3개만 쓰면 양이온이 donor 축을 중심으로 돌 수 있었음: c4 config에서 경로 이탈 최대 2.94 → 0.52 Å (2026-11-05, 우리 계산).
  - 양이온이 의도한 점에서 1.0 Å 넘게 벗어난 frame(예: 기준 평면에 대한 거울상)은 버림. params `path_dev_A`, `ref3`.
  - 양이온이 어떤 원자에 공유결합 반지름 합의 0.75배보다 가까워진 frame은 버림.
  - 양이온–donor constraint 거리 자체가 공유결합 반지름 합의 0.75배보다 짧은 점은 xTB를 돌리지 않고 건너뜀.
  - `exchange: true`면 첫 배위권 donor 하나를 양이온에서 `exchange_extend_A`(기본 3.0 Å)만큼 멀어지게 하는 relaxed scan도 만듦 (ligand exchange의 이탈 쪽). 양이온이 다른 원자와 겹치는 frame은 버림. params `kind: exchange_out`. cluster당 `exchange_max_per_cluster`개.

## Charge와 multiplicity

- charge = 구성 성분 형식 전하의 합 + `delta_q`.
- 최저 multiplicity는 전자 수가 짝수면 1, 홀수면 2. 여기에 `mult_offsets`를 더함.
- `bs: true`인 sampler는 전자 수가 짝수일 때 singlet을 BS(broken-symmetry) singlet(`guess=mix`)으로만 만듦. homolytic 절단 경로를 위한 것임.
- xTB 경로는 charge마다 따로 계산함 (uhf = 전자 수 % 2). 같은 경로를 여러 multiplicity로 계산함.
- 환원(Δq < 0)은 전체 전하가 −1 이상이면서 Li/Na가 있거나 중성 radical(예: CH3O)이 있는 계에만 만듦. 양이온 없는 용매·음이온의 추가 전자나 dianion은 기체상에서 결합하지 않아 diffuse basis 꼬리에 갇힌 가짜 상태가 됨 (Li + 2 PF6의 q = −2도 Li⁰ + 2 PF6⁻라 제외). GFN2-xTB는 이를 판별하지 못함: H2O⁻ vertical detachment energy +3.8 eV, DME⁻ +5.7 eV (2026-11-05, 우리 계산).
- Δq ≤ −2는 Li/Na를 가진 성분 수가 |Δq| 이상인 계에만 만듦 (전자 하나당 Li/Na 하나).
- 원자 겹침 판정은 f × (공유결합 반지름 합)과 0.6 Å 중 큰 값을 씀 (H–H는 0.7배가 0.43 Å라 하한을 둠).

## Gaussian 입력

```
%nprocshared=16
%mem=32GB
#p {R|U}wB97XD/gen nosymm int=ultrafine scf=(xqc,maxcycle=300) pop=(npa,hirshfeld) iop(6/80=1) force | pop=(npa,hirshfeld,always) iop(6/80=1) opt=(maxcycles=200) [guess=mix]

<id>

<charge> <mult>
<좌표>

<원소별 def2-TZVPPD block>
```

- `R`은 closed-shell singlet에 씀. multiplicity > 1이거나 BS singlet이면 `U`.
- `nosymm`을 써서 force가 입력 좌표계로 나옴.
- `pop=(npa,hirshfeld)`는 NBO natural charge와 Hirshfeld charge·spin, CM5 charge를, `iop(6/80=1)`은 Loewdin population을 출력하게 함. opt job은 `pop=(npa,hirshfeld,always)`로 매 step population을 출력함 (`g16.pop_always: false`면 `pop=(npa,hirshfeld)`).
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
