/**
 * Online typing via RIVM MPF Typing Tools (Playwright automation).
 * Supports both Norovirus and Rotavirus A.
 */

const fs = require('fs');
const { chromium } = require('playwright');

const TOOLS = {
  noro: {
    url: 'https://mpf.rivm.nl/mpf/typingtool/norovirus/',
    label: 'Norovirus',
  },
  rota: {
    url: 'https://mpf.rivm.nl/mpf/typingtool/rotavirusa/',
    label: 'Rotavirus A',
  },
};

// Simple CSV parser
function parseCSV(csvText) {
  const lines = csvText.trim().split('\n');
  if (lines.length < 2) return [];
  
  const parseLine = (line) => {
    const result = [];
    let current = '';
    let inQuotes = false;
    for (let i = 0; i < line.length; i++) {
      const char = line[i];
      if (char === '"') {
        inQuotes = !inQuotes;
      } else if (char === ',' && !inQuotes) {
        result.push(current.trim());
        current = '';
      } else {
        current += char;
      }
    }
    result.push(current.trim());
    return result.map(v => v.replace(/^"|"$/g, ''));
  };

  const headers = parseLine(lines[0]);
  const rows = [];
  for (let i = 1; i < lines.length; i++) {
    if (!lines[i].trim()) continue;
    const values = parseLine(lines[i]);
    const row = {};
    headers.forEach((h, index) => {
      row[h] = values[index] || '';
    });
    rows.push(row);
  }
  return rows;
}

/**
 * Submit FASTA to a single RIVM typing tool and return results.
 */
async function submitToRIVM(virus, fastaContent, opts = {}) {
  const tool = TOOLS[virus];
  if (!tool) throw new Error(`Unknown virus: ${virus}`);

  const { proxy, timeout = 360000 } = opts;

  const launchOpts = { headless: true };
  if (proxy) launchOpts.proxy = { server: proxy };

  const browser = await chromium.launch(launchOpts);
  const context = await browser.newContext({
    userAgent: 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/131.0.0.0 Safari/537.36',
  });
  const page = await context.newPage();

  try {
    // 1. Navigate
    await page.goto(tool.url, { waitUntil: 'domcontentloaded', timeout: 30000 });
    await page.waitForSelector('textarea', { timeout: 15000 });

    // 2. Paste FASTA
    const textarea = page.locator('textarea');
    await textarea.fill(fastaContent);

    // 3. Click Start
    await page.evaluate(() => {
      const btn = document.querySelector('button[id^="button-run"]');
      if (btn) btn.click();
    });

    // 4. Handle Warning dialog
    try {
      const yesBtn = page.locator('button:has-text("Yes")');
      await yesBtn.waitFor({ state: 'visible', timeout: 5000 });
      await yesBtn.click();
    } catch {
      // No warning dialog — proceed
    }

    // 5. Wait for CSV download link to appear (indicates completion)
    const csvLink = page.locator('a:has-text("CSV")');
    await csvLink.waitFor({ state: 'visible', timeout });

    // Extract Job ID if possible
    let jobId = '';
    try {
      jobId = await page.evaluate(() => {
        const links = document.querySelectorAll('a');
        for (const a of links) {
          const m = a.textContent.match(/Monitor job \[(\d+)\]/);
          if (m) return m[1];
        }
        return '';
      });
    } catch {
      // Could not detect job ID
    }

    // Get CSV URL
    const csvUrl = await csvLink.getAttribute('href') || '';

    // 6. Download CSV
    const [ download ] = await Promise.all([
      page.waitForEvent('download'),
      csvLink.click(),
    ]);

    const downloadPath = await download.path();
    const csvContent = fs.readFileSync(downloadPath, 'utf8');

    // Parse CSV to results
    const results = parseCSV(csvContent);

    return { virus, jobId, results, csvUrl, csvContent };
  } finally {
    await browser.close();
  }
}

/**
 * Run both norovirus and rotavirus online typing in parallel.
 * Returns { noro: {...}, rota: {...} }
 */
async function runOnlineTyping(fastaContent, opts = {}) {
  const [noro, rota] = await Promise.allSettled([
    submitToRIVM('noro', fastaContent, opts),
    submitToRIVM('rota', fastaContent, opts),
  ]);

  return {
    noro: noro.status === 'fulfilled' ? noro.value : { error: noro.reason?.message },
    rota: rota.status === 'fulfilled' ? rota.value : { error: rota.reason?.message },
  };
}

module.exports = { submitToRIVM, runOnlineTyping };
