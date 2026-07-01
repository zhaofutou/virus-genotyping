# Virus Genotyping — Online + Local

Unified virus subtyping tool combining **online (RIVM MPF)** and **local (reference-based)** methods for Norovirus and Rotavirus A.

## Features

- **Online typing**: Automated submission to RIVM typing tools via Playwright
  - Norovirus: https://mpf.rivm.nl/mpf/typingtool/norovirus/
  - Rotavirus A: https://mpf.rivm.nl/mpf/typingtool/rotavirusa/
- **Local typing**: Reference-based alignment using k-mer indexing + ungapped overlap
- **Parallel execution**: Online and local typing run simultaneously
- **Monthly convention support**: Files named `May05V4`, `Jun12GII`, etc. get a 1-72 index summary
- **Side-by-side comparison**: Online vs local results displayed together
- **ZIP upload**: Upload individual FASTA files or a ZIP archive

## Quick Start

```bash
git clone <repo-url> && cd virus-genotyping
npm install
npx playwright install chromium
node server.js
```

Open http://localhost:3000 and upload your FASTA files.

## API

```bash
# Upload files
curl -X POST http://localhost:3000/api/upload \
  -F "files=@sample1.fasta" -F "files=@sample2.fasta"

# Check status (SSE)
curl http://localhost:3000/api/status/<jobId>

# Download results ZIP
curl http://localhost:3000/api/download/<jobId> -o results.zip
```

## Naming Convention

For 1-72 summary tables, name files following the monthly pattern:

| Pattern | Meaning |
|---------|---------|
| `May05V4` | Index 5, Rotavirus VP4 |
| `May05V7` | Index 5, Rotavirus VP7 |
| `Jun12GI` | Index 12, Norovirus GI |
| `Jun12GII` | Index 12, Norovirus GII |

## Output Files

| File | Description |
|------|-------------|
| `online_noro_result.csv` | RIVM Norovirus results |
| `online_rota_result.csv` | RIVM Rotavirus results |
| `local_noro_result.csv` | Local Norovirus results |
| `local_rota_result.csv` | Local Rotavirus results |
| `online_final_1_72.csv` | Online 1-72 summary (if naming convention used) |
| `local_final_1_72.csv` | Local 1-72 summary (if naming convention used) |

## Tech Stack

- Node.js + Express
- Playwright (headless Chromium for RIVM automation)
- Worker threads (parallel local typing)
- AdmZip (ZIP handling)
