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

- key 이름은 `6_build_charge-aware-mlip`과 같은 규칙을 따름.
- spin contamination(`s2_dev`)이 큰 frame도 값만 저장하고 걸러내지 않음.
- NBO spin은 NBO 출력의 alpha와 beta natural population 차이로 계산함.

## 확인 결과 (Hangang G16 c02, job 1212786, 2026-11-04, test_decks/ t1–t3)

- Loewdin: header는 `Lowdin Atomic Charges:`. open-shell에서도 charge만 출력하고 spin은 출력하지 않음 → open-shell `lowdin_spins`는 NaN.
- NPA: total / Alpha / Beta block 형식 확인. NBO spin = α − β natural population. 합은 doublet 1.000, BS singlet 0.000.
- opt: 기본은 population이 첫 step과 수렴 후에만 출력됨 (t1). `pop=(npa,always)`면 매 step 출력됨: t4(t1과 같은 deck, Hangang job 1217563, 2026-11-05)에서 8개 frame 모두 Mulliken, Loewdin, NBO charge와 NBO spin이 있음. 시간 199 s (t1 157 s, 다른 node라 참고값).
- restricted closed-shell은 세 방식 모두 spin이 0.
- 주의: def2-TZVPPD의 diffuse 함수 때문에 Mulliken/Loewdin 값이 basis에 따라 크게 바뀜.
  - 예: Li⁰·EC의 Mulliken spin이 Li 1.499 (NBO 0.926).
  - 예: Li+–EC의 Loewdin carbonyl O charge가 +0.09.
