let omicsqcSelectionNonce = null;
let omicsqcSelectionCounter = 0;
document.addEventListener("change", function (event) {
  const input = event.target;
  if (input.id !== "assay" || !input.files || !input.files.length) return;

  const oversize = input.files[0].size > 500 * 1024 * 1024;
  omicsqcSelectionNonce = String(++omicsqcSelectionCounter);
  const run = document.getElementById("run");
  if (run) run.disabled = true;
  Shiny.setInputValue("assay_selection", { oversize: oversize, nonce: omicsqcSelectionNonce }, { priority: "event" });

  if (oversize) {
    event.stopImmediatePropagation();
    input.value = "";
    Shiny.setInputValue("assay", null, { priority: "event" });
  }
}, true);

Shiny.addCustomMessageHandler("omicsqc-run-button", function (state) {
  const run = document.getElementById("run");
  if (run) run.disabled = !(state.enabled && state.nonce === omicsqcSelectionNonce);
});

Shiny.addCustomMessageHandler("omicsqc-redirect", function (url) {
  const target = document.documentElement.classList.contains("cvdlink-embedded") ? window.top : window;
  target.location.replace(url);
});

Shiny.addCustomMessageHandler("omicsqc-disable-upload", function (disabled) {
  const input = document.getElementById("assay");
  if (input) input.disabled = disabled;
});

let omicsqcCookieName = null;
Shiny.addCustomMessageHandler("omicsqc-session-cookie", function (cookie) {
  omicsqcCookieName = cookie.name;
  document.cookie = cookie.name + "=" + cookie.value + "; Path=/; SameSite=Strict";
});

window.addEventListener("pagehide", function () {
  if (omicsqcCookieName) {
    document.cookie = omicsqcCookieName + "=; Path=/; Max-Age=0; SameSite=Strict";
  }
});
