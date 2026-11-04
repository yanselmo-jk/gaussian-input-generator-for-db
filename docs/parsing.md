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
5. opt job은 최적화의 모든 step이 frame이 됨. Gaussian은 charge를 첫 step과 수렴 후에만 출력하므로 중간 step의 charge는 NaN임.

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

## 미검증 사항

`test_decks/`의 Gaussian 입력 3개를 실행한 log로 확인한 뒤 parser를 확정할 예정임.

| test deck | 확인할 것 |
|---|---|
| `t1_li_ec_q1_m1_opt.gjf` | closed-shell opt. Loewdin/NPA charge 출력 형식, opt 중간 step의 population 출력 여부 |
| `t2_li_ec_q0_m2_sp.gjf` | 환원 doublet. open-shell NPA alpha/beta block 형식, Loewdin spin 출력 형식 |
| `t3_ec_co_stretch_bs_sp.gjf` | 고리 C–O를 늘인 BS singlet (`guess=mix`) |

1. `iop(6/80=1)`이 Loewdin charge/spin을 출력하는지와 그 header 문구. 지금은 `Lowdin Atomic Charges` 형태로 가정함.
2. open-shell NPA의 alpha/beta block 형식.
3. opt 중간 step에 population이 출력되는지.
