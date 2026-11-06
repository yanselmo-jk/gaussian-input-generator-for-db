# Production 작업 순서 (두 machine)

생성 machine(A)에서 입력을 만들고, 계산 machine(B)에서 Gaussian16을 돌리고, A에서 결과를 학습 data와 DB로 모으는 순서.
명령은 모두 A의 프로젝트 폴더(`rmd generate`를 실행한 폴더)에서 실행함.

## 0. 준비 (A)

```bash
conda activate rmlipdata
rmd species
```

- config는 `configs/example.yaml`을 복사해 고침. 화학종(`species`), charge 변화(`charge.delta_q`), sampler별 개수를 정함.
- `out:`이 이 계산 묶음 전체의 폴더임 (예: `runs/prod1`). 화학종을 늘려 가며 여러 번 생성할 때도 같은 `out`을 씀.

## 1. 생성 (A)

```bash
rmd generate runs/prod1.yaml --dry
rmd generate runs/prod1.yaml
```

- `--dry`: xTB 없이 sampler별 생성 단위 수만 출력함. 규모를 먼저 확인함.
- 결과: `runs/prod1/gjf/<sp|opt>/<sampler>/<id>.gjf`, `runs/prod1/registry.sqlite`(DB), `runs/prod1/manifest.csv`(아직 보내지 않은 deck 목록).
- 중간에 멈추면 같은 명령을 다시 실행함. 이미 쓴 deck은 다시 쓰지 않음.
- 계산 자원을 바꾸려면 보내기 전에:

```bash
rmd resources --registry runs/prod1/registry.sqlite --nproc 32 --mem 64GB
```

## 2. 묶어서 보내기 (A → B)

```bash
rmd pack bundle_01 --registry runs/prod1/registry.sqlite
scp bundle_01.tar.gz USER@HOST_B:/path/to/work/
```

- `bundle_01/<id>/in.gjf`, `bundle_01/jobs.csv`(id, nproc, mem), `bundle_01.tar.gz`가 생김.
- 묶은 deck은 DB에서 `submitted`가 됨. 다음 pack에는 그 뒤 새로 생긴 deck만 들어감.
- sampler별로 나눠 묶으려면 `--sampler cluster pair_scan`.
- 잘못 묶었으면 `rmd submit bundle_01 --registry runs/prod1/registry.sqlite --undo`.

## 3. 계산 (B)

```bash
tar xzf bundle_01.tar.gz
cd bundle_01
for d in */; do sbatch --chdir="$d" ~/job_g16.sh; done
```

- `~/job_g16.sh` 안에서 `g16 in.gjf > in.log`를 실행함. 폴더 하나가 job 하나.
- nproc, mem은 각 `in.gjf` 첫 줄과 `jobs.csv`에 있음.

## 4. 결과 가져오기와 parse (B → A)

```bash
rsync -a --include='*/' --include='in.log' --exclude='*' USER@HOST_B:/path/to/work/bundle_01/ back/bundle_01/
rmd parse back/ -o data.extxyz --registry runs/prod1/registry.sqlite
```

- `back/` 아래 모든 묶음의 `in.log`를 한 번에 읽음. id는 log 안의 deck title에서 읽으므로 폴더·파일 이름과 무관함.
- 계산 중인 job이 섞여 있어도 됨: 종료 줄이 없는 log는 건너뛰고 상태를 바꾸지 않음. 나중에 다시 가져와 parse하면 끝난 job이 추가됨.
- `data.extxyz`는 매번 새로 씀. 항상 `back/` 전체를 parse함 (같은 경로의 frame을 함께 비교해야 `path_flag`가 맞음).
- 상태: 정상 종료 → `parsed`, `Error termination` → `failed`.

## 5. 진행 상황과 재계산 (A)

```bash
rmd status --registry runs/prod1/registry.sqlite
rmd list --registry runs/prod1/registry.sqlite --status failed
rmd list --registry runs/prod1/registry.sqlite --status submitted
```

- `failed`: Gaussian 오류. 원인(log)을 보고 다시 돌릴 deck은 B의 같은 폴더에서 다시 제출하고 A에서 `rmd submit <id> --registry ...`로 `submitted`로 되돌림.
- 오래 `submitted`로 남은 deck: 계산 중이거나 walltime으로 끊긴 job. B에서 확인 후 다시 제출함.

## 6. 화학종을 늘려 다시 생성 (A)

```bash
rmd generate runs/prod1_more.yaml
rmd pack bundle_02 --registry runs/prod1/registry.sqlite
```

- `prod1_more.yaml`은 `species`만 늘리고 `out`, `seed`, sampler 설정은 그대로 둠.
- 이미 만든 조합은 같은 구조(같은 id)가 나와 다시 쓰지 않음. bundle_02에는 새 deck만 들어감.
- sampler 설정(`n_points`, `n_distort` 등)이나 seed를 바꾸면 같은 조합이라도 새 구조가 생김.
- `max_combos`, `max_pairs` 상한 안에서 조합을 다시 뽑으므로, 이미 한 조합이 자리를 차지해 새 조합이 기대보다 적을 수 있음. config에 `cluster: {skip_done: true}`를 두면 이미 한 조합(현재 delta_q 전부에서 기록이 있는 것)을 후보에서 빼고 상한을 새 조합에 씀. delta_q를 새로 추가했다면 기존 조합도 그 charge가 없으므로 다시 뽑힘.

## 7. 학습 data (A)

- `data.extxyz`: frame마다 energy(eV), forces(eV/Å), charge, spin(= multiplicity), Mulliken/Löwdin/NPA/Hirshfeld charge·spin, CM5 charge, dipole_eA, quadrupole_eA2, homo_eV, s2, sampler, composition, group, path_flag 등 (전체 key: docs/parsing.md).
- 이미 거른 frame: 원자 force > 50 eV/Å, ⟨S²⟩ − S(S+1) 초과(평형 sampler의 BS 아닌 frame 0.25, 그 밖 1.1), HOMO > 0 eV. `--keep-all`이면 거르지 않음.
- 표시만 한 것: `path_flag`(경로 에너지·spin 불연속), `q_sum_ok`, `spin_sum_ok`, `converged`. 학습 전에 직접 고름.
- train/test는 `group`(= 조성)으로 나누면 같은 cluster·경로가 양쪽에 들어가지 않음.

## 8. DB 조회 (A)

```bash
rmd export --registry runs/prod1/registry.sqlite -o registry.csv
sqlite3 runs/prod1/registry.sqlite "SELECT sampler, status, COUNT(*) FROM structures GROUP BY 1, 2"
```

- 표 구조와 예시 질의: docs/registry.md.

## 주의

- 한 계산 묶음에는 하나의 `out` 폴더를 씀. 다른 폴더에 생성하면 DB가 따로라 같은 구조도 다시 만듦.
- `rmd generate`가 도는 동안 같은 DB에 parse·pack을 하지 않음 (최대 10분 기다린 뒤 멈출 수 있음).
- 결과를 Windows를 거쳐 옮겨도 됨 (CRLF log도 읽음).
