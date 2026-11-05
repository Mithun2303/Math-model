// Reads JSON [{tex, display}] from argv[2], writes JSON [html] to argv[3]
const katex = require('/tmp/pdftool/node_modules/katex');
const fs = require('fs');
const items = JSON.parse(fs.readFileSync(process.argv[2], 'utf8'));
const out = items.map(it => {
  try { return {ok: true, html: katex.renderToString(it.tex, {displayMode: it.display, throwOnError: true, strict: 'ignore', trust: false})}; }
  catch (e) { return {ok: false, err: String(e.message).slice(0, 300), tex: it.tex}; }
});
fs.writeFileSync(process.argv[3], JSON.stringify(out));
