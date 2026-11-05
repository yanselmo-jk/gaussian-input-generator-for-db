# Registry (SQLite) 사용법

## SQLite란

- SQLite는 표(table)를 파일 하나에 담는 database임.
- 서버가 필요 없고, `registry.sqlite` 파일을 복사하면 그대로 옮겨짐.
- 이 파일에는 `structures`라는 표가 하나 있음. 한 행이 Gaussian 입력 하나에 해당함.

| column | 의미 |
|---|---|
| `id` | 구조 id (= `.gjf`/`.log` 파일명) |
| `parent_id` | vertical deck(같은 좌표, charge ±1)의 원래 deck id (최저 multiplicity deck). 그 밖에는 비어 있음 |
| `sampler` | 만든 sampler |
| `composition` | 구성 (예: `Li1 PF6_1 EC2`) |
| `formula`, `n_atoms` | 화학식, 원자 수 |
| `charge`, `mult` | 전하, multiplicity |
| `bs` | BS singlet 여부 (0/1) |
| `job` | `sp` 또는 `opt` |
| `params` | sampler별 상세 정보. JSON 문자열 (scan 좌표, 경로 index, template, `is_peak` 등). xTB SCC가 수렴하지 않아 `--etemp 1000`으로 다시 돌린 구조에는 `xtb_etemp_K: 1000`이 붙음 |
| `gjf_path`, `log_path` | 입력 파일과 log 파일 경로 |
| `status` | `generated` (입력만 있음) → `submitted` (`rmd submit`으로 표시함) → `parsed` (log 읽음) 또는 `failed` (비정상 종료) |
| `created` | 생성 시각 |

## 보는 방법

1. **CSV로 내보내기**: 가장 쉬운 방법. Excel 등으로 열면 됨.
   ```bash
   rmd export --registry runs/example/registry.sqlite -o registry.csv
   ```
2. **요약 보기**:
   ```bash
   rmd status --registry runs/example/registry.sqlite
   ```
3. **터미널에서 직접 질의(query)**:
   ```bash
   sqlite3 runs/example/registry.sqlite
   sqlite> .headers on
   sqlite> .mode column
   sqlite> SELECT sampler, status, COUNT(*) FROM structures GROUP BY sampler, status;
   sqlite> .quit
   ```
4. **GUI**: [DB Browser for SQLite](https://sqlitebrowser.org/)로 파일을 열면 표를 마우스로 볼 수 있음.

## 자주 쓰는 query

```sql
-- 계산 실패한 구조
SELECT id, log_path FROM structures WHERE status='failed';

-- 구성, charge, multiplicity별 개수
SELECT composition, charge, mult, COUNT(*) FROM structures GROUP BY composition, charge, mult;

-- BS singlet인 reaction 구조 5개
SELECT id, composition, charge FROM structures WHERE sampler='reaction' AND bs=1 LIMIT 5;

-- params(JSON) 안의 값 꺼내기: reaction의 pseudo-TS 최고점만
SELECT id, json_extract(params, '$.template') AS template FROM structures
WHERE sampler='reaction' AND json_extract(params, '$.is_peak')=1;
```

## Python에서 extxyz와 합치기 (raw data DB를 만들 때)

extxyz frame의 `info["id"]`가 registry의 `id`와 같음. 이것으로 연결함.

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
