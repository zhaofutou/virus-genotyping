/**
 * Online typing via RIVM MPF Typing Tools (Playwright automation).
 * Supports both Norovirus and Rotavirus A.
 */

const { chromium } = require('playwright');

const TOOLS = {
  noro: {
    url: 'https://mpf.rivm.nl/mpf/typingtool/norovirus/',
    label: 'Norovirus',
    pollInterval: 6000,
    maxPoll: 60,
  },
  rota: {
    url: 'https://mpf.rivm.nl/mpf/typingtool/rotavirusa/',
    label: 'Rotavirus A',
    pollInterval: 6000,
    maxPoll: 60,
  },
};

/**
 * Submit FASTA to a single RIVM typing tool and return results.
 */
async function submitToRIVM(virus, fastaContent, opts = {}) {
  const tool = TOOLS[virus];
  if (!tool) throw new Error(`Unknown virus: ${virus}`);

  const { proxy, timeout = 360000 } = opts;
  const maxPoll = Math.ceil(timeout / tool.pollInterval);

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
      await yesBtn.waitFor({ state: 'visible', timeout: 10000 });
      await yesBtn.click();
    } catch {
      // No warning dialog — proceed
    }

    // 5. Extract Job ID
    let jobId = '';
    try {
      await page.waitForFunction(
        () => {
          const links = document.querySelectorAll('a');
          return Array.from(links).some(a => /Monitor job \[\d+\]/.test(a.textContent));
        },
        { timeout: 15000 }
      );
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

    // 6. Poll for results
    let results = [];
    let csvUrl = '';

    for (let attempt = 0; attempt < maxPoll; attempt++) {
      await page.waitForTimeout(tool.pollInterval);

      const tableData = await page.evaluate(() => {
        const table = document.querySelector('table');
        if (!table) return null;

        const tbodyRows = table.querySelectorAll('tbody tr');
        if (tbodyRows.length === 0) return null;

        const headers = Array.from(table.querySelectorAll('thead th')).map(th => ({
          text: th.textContent.trim(),
          colspan: th.colSpan,
        }));

        // Expand headers: colspan=2 means "type + subtype"
        const fixedHeaders = [];
        for (const h of headers) {
          if (h.colspan === 2) {
            fixedHeaders.push(h.text + ' type');
            fixedHeaders.push(h.text + ' subtype');
          } else {
            fixedHeaders.push(h.text);
          }
        }

        const rows = [];
        for (const tr of tbodyRows) {
          const cells = Array.from(tr.children);
          const rowData = {};
          cells.forEach((cell, i) => {
            rowData[fixedHeaders[i] || `col_${i}`] = cell.textContent.trim();
          });
          rows.push(rowData);
        }

        let csv = '';
        const links = document.querySelectorAll('a');
        for (const a of links) {
          if (a.textContent.includes('CSV')) {
            csv = a.href;
            break;
          }
        }

        return { headers: fixedHeaders, rows, csvUrl: csv };
      });

      if (tableData && tableData.rows.length > 0) {
        results = tableData.rows;
        csvUrl = tableData.csvUrl;
        break;
      }
    }

    if (results.length === 0) {
      throw new Error(`${tool.label}: No results after ${maxPoll} polling attempts (${(maxPoll * tool.pollInterval) / 1000}s).`);
    }

    // 7. Download CSV
    let csvContent = '';
    if (csvUrl) {
      try {
        csvContent = await page.evaluate(async (url) => {
          const resp = await fetch(url);
          return await resp.text();
        }, csvUrl);
      } catch {
        // CSV download failed
      }
    }

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
