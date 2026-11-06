# Log parser (`rmd parse`)

```bash
rmd parse logs/ -o data.extxyz --registry runs/example/registry.sqlite
```

- PATH에는 log 파일이나 폴더를 줌. 폴더면 안의 `*.log`, `*.out`을 모두 읽음.
- `--registry`를 주면 각 log의 status를 바꿈. log 파일명(= id)으로 행을 찾음.
  - 모든 job step이 정상 종료되고 frame이 있으면 `parsed`, 아니면 `failed`.
- 순서: extxyz를 임시 파일에 다 쓰고 교체한 다음 registry 변경을 확정함. 그래서 저장에 실패하면 registry가 바뀌지 않음.

## Frame을 나누는 규칙

1. log를 Link1 job step 단위로 나눔.
2. 각 step 안에서 `Input orientation:`이 나올 때마다 새 frame을 시작함.
3. `SCF Done`과 force 표가 둘 다 있는 부분만 frame이 됨.
4. 최적화가 끝난 뒤(`Stationary point found`) 다시 출력되는 population은 같은 좌표의 마지막 frame에 합침.
5. opt job은 최적화의 모든 step이 frame이 됨. 생성기가 쓰는 opt deck은 기본값(`g16.pop_always: true`)에서 `pop=(npa,always)`라서 매 step의 charge와 spin이 있음. `always` 없이 돌린 opt log는 첫 step과 수렴 후에만 population이 있어 중간 step은 NaN임.

## extxyz key

| 구분 | key | 단위 / 의미 |
|---|---|---|
| 좌표 | `positions` | Å |
| energy | `energy` | eV (SCF energy × 27.211386245988) |
| force | `forces` | eV/Å (Hartree/Bohr × 27.211386245988 / 0.529177210903) |
| charge | `mulliken_charges`, `lowdin_charges`, `nbo_charges` | e. 해당 frame에 출력이 없으면 NaN |
| spin | `mulliken_spins`, `lowdin_spins`, `nbo_spins` | α−β spin density. restricted closed-shell은 0 |
| info | `id`, `step`, `frame` | registry id, Link1 step, frame 번호 |
| info | `charge`, `spin` | 전체 전하, multiplicity (`spin` = multiplicity) |
| info | `job`, `converged`, `normal_termination` | job 종류, 수렴 여부(opt 마지막 frame), 정상 종료 여부 |
| info | `s2`, `s2_dev` | ⟨S²⟩, ⟨S²⟩ − S(S+1) |
| info | `q_sum_ok`, `spin_sum_ok` | Mulliken charge 합 = 전체 전하 (±1e-3), spin 합 = mult − 1 (±0.05) |
| charge | `hirshfeld_charges`, `cm5_charges` | e. `pop=hirshfeld` 출력. 출력이 없으면 NaN |
| spin | `hirshfeld_spins` | Hirshfeld spin density. closed shell은 0 |
| info | `dipole` | e·Å, (x, y, z). Gaussian "Dipole moment (field-independent basis)" (Debye) × 0.20819434 |
| info | `quadrupole` | e·Å², traceless (XX, YY, ZZ, XY, XZ, YZ). Gaussian "Traceless Quadrupole moment" (Debye·Å) × 0.20819434 |
| info | `homo_eV` | 가장 높은 점유 궤도 에너지 (α와 β 중 큰 값), eV |
| info | `bs` | route에 `guess=mix`가 있으면 true (BS singlet) |
| info | `path`, `path_pos`, `rattled`, `path_flag` | 경로 key, 경로 위치, 흔든 frame 여부, 연속성 표시 (`--registry`를 줄 때, 경로 sp frame만) |
| info | `sampler`, `composition`, `group` | `--registry`를 줄 때만. `group` = 조성. train/test를 group으로 나누면 같은 cluster·경로가 양쪽에 들어가지 않음 |

- key 이름은 `6_build_charge-aware-mlip`과 같은 규칙을 따름.
- `rmd parse`는 다음 frame을 extxyz에 쓰지 않고 개수만 출력함 (`--keep-all`이면 거르지 않음). registry status는 필터와 상관없이 log 단위로 정해짐.
  - 원자 하나의 force가 50 eV/Å 초과 (OMol25 Appendix A.1).
  - `s2_dev`가 BS singlet은 1.1, 그 밖의 open-shell은 0.25 초과 (Sc–Zn이 있으면 0.5와 비교해 작은 값). 1.1은 OMol25 기준을 이상값과의 차이로 읽은 값(절대값이면 triplet이 모두 빠짐)이고 BS singlet의 완전한 spin 분리(차이 약 1)를 허용함. 0.25는 우리 선택임 (예: doublet인데 ⟨S²⟩ 1.8이면 오염된 상태). Li/Na는 metal 기준에 넣지 않음.
  - `homo_eV` > 0: 추가 전자가 결합하지 않고 diffuse basis 꼬리에 있음.
- NBO spin은 NBO 출력의 alpha와 beta natural population 차이로 계산함.

## 경로 연속성 표시 (`path_flag`)

`--registry`를 주면 같은 경로의 sp frame을 경로 순서대로 비교해 `path_flag`를 붙임. 같은 경로는 sampler, 조성, charge, multiplicity, BS, 경로 key가 같은 frame임. 경로 key는 bond_scan의 bond, reaction의 template과 원자 쌍, hop의 양이온과 donor임. frame은 버리지 않음.

- `ok`: 이상 없음.
- `energy`: 흔들지 않은 이웃 frame 사이 원자 최대 이동이 0.5 Å 미만인데, 에너지 변화가 force로 예상한 값 −½(F_a + F_b)·(R_b − R_a)와 2 eV 넘게 다름. SCF가 다른 해에 수렴했을 가능성이 있음.
- `spin`: 원자 하나의 NPA spin이 이웃 frame 사이에서 0.5 e 넘게 바뀜. 홀전자가 옮겨 감. 실제 반응(H 이동, 결합 절단)에서도 생길 수 있음.
- 기준값은 우리 선택임. c4 test config의 xTB 경로로 시험한 결과(2026-11-06, 우리 계산)를 근거로 정함. 이동이 크거나 흔든 frame이 낀 구간은 예상값 오차가 커서(중앙값 0.4–2.6 eV) 에너지 검사에서 뺌.
- 단원자, vertical deck, 최고점 rattle 복사본은 경로에 넣지 않음.

## 확인 결과 (Hangang G16 c02, job 1212786, 2026-11-04, test_decks/ t1–t3)

- Loewdin: header는 `Lowdin Atomic Charges:`. open-shell에서도 charge만 출력하고 spin은 출력하지 않음 → open-shell `lowdin_spins`는 NaN.
- NPA: total / Alpha / Beta block 형식 확인. NBO spin = α − β natural population. 합은 doublet 1.000, BS singlet 0.000.
- opt: 기본은 population이 첫 step과 수렴 후에만 출력됨 (t1). `pop=(npa,always)`면 매 step 출력됨: t4(t1과 같은 deck, Hangang job 1217563, 2026-11-05)에서 8개 frame 모두 Mulliken, Loewdin, NBO charge와 NBO spin이 있음. 시간 199 s (t1 157 s, 다른 node라 참고값).
- restricted closed-shell은 세 방식 모두 spin이 0.
- 주의: def2-TZVPPD의 diffuse 함수 때문에 Mulliken/Loewdin 값이 basis에 따라 크게 바뀜.
  - 예: Li⁰·EC의 Mulliken spin이 Li 1.499 (NBO 0.926).
  - 예: Li+–EC의 Loewdin carbonyl O charge가 +0.09.

## 확인 결과 (Hangang G16 c02, job 1224122, 2026-11-06)

- `pop=(npa,hirshfeld)` 출력 형식: EC⁺ doublet sp, Li⁺–EC opt (`pop=(npa,hirshfeld,always)`에서 매 step 출력). fixture `tests/fixtures/h1_ec_cation_hirsh.log`.
- 같은 job에서 `stable=opt` + Link1 force도 시험함: EC⁺ doublet은 stable. 안정성 검사 92 s(16 core) 대 sp + force 45 s. deck에는 넣지 않음.
