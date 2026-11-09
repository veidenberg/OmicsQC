const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const test = require("node:test");
const vm = require("node:vm");

test("the Run button follows the current upload and exact size limit", () => {
  const handlers = {};
  const inputs = [];
  const run = { disabled: false };
  const assay = { id: "assay", files: [], value: "selected" };
  let onChange;
  const context = {
    document: {
      addEventListener(name, callback) {
        if (name === "change") onChange = callback;
      },
      getElementById(id) {
        return id === "run" ? run : assay;
      }
    },
    window: { addEventListener() {} },
    Shiny: {
      addCustomMessageHandler(name, callback) { handlers[name] = callback; },
      setInputValue(name, value) { inputs.push({ name, value }); }
    }
  };
  vm.runInNewContext(fs.readFileSync(path.join(__dirname, "../www/upload-limit.js"), "utf8"), context);

  function select(size) {
    assay.files = [{ size }];
    let stopped = false;
    onChange({ target: assay, stopImmediatePropagation() { stopped = true; } });
    const selection = inputs.filter(input => input.name === "assay_selection").at(-1).value;
    return { selection, stopped };
  }

  const first = select(500 * 1024 * 1024);
  assert.equal(first.selection.oversize, false);
  assert.equal(first.stopped, false);
  assert.equal(run.disabled, true);
  handlers["omicsqc-run-button"]({ enabled: true, nonce: first.selection.nonce });
  assert.equal(run.disabled, false);

  const next = select(1);
  handlers["omicsqc-run-button"]({ enabled: true, nonce: first.selection.nonce });
  assert.equal(run.disabled, true);
  handlers["omicsqc-run-button"]({ enabled: true, nonce: next.selection.nonce });
  assert.equal(run.disabled, false);

  const tooLarge = select(500 * 1024 * 1024 + 1);
  assert.equal(tooLarge.selection.oversize, true);
  assert.equal(tooLarge.stopped, true);
  assert.equal(run.disabled, true);
  assert.equal(assay.value, "");
  assert.equal(inputs.at(-1).name, "assay");
  assert.equal(inputs.at(-1).value, null);
});

test("authentication redirects leave the dashboard frame when embedded", () => {
  const handlers = {};
  const redirects = [];
  let embedded = false;
  const context = {
    document: {
      documentElement: { classList: { contains: () => embedded } },
      addEventListener() {}
    },
    window: {
      addEventListener() {},
      location: { replace: url => redirects.push(["window", url]) },
      top: { location: { replace: url => redirects.push(["top", url]) } }
    },
    Shiny: { addCustomMessageHandler(name, callback) { handlers[name] = callback; } }
  };
  vm.runInNewContext(fs.readFileSync(path.join(__dirname, "../www/upload-limit.js"), "utf8"), context);
  handlers["omicsqc-redirect"]("/login");
  embedded = true;
  handlers["omicsqc-redirect"]("/login");
  assert.deepEqual(redirects, [["window", "/login"], ["top", "/login"]]);
});
