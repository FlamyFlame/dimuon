#!/usr/bin/env python3
"""No-muon-pT-filter signal-acceptance sample (signal_acceptance_mass_window.md D5). Copy of
   evgen_high_pT_additional/grid/chain.py (tested 2026-09); only names/paths differ.
   Drive the grid chain per pT-hat slice:
   evgen (prun, already submitted) -> Derivation_tf TRUTH0 (pathena, Athena 25.0.72)
   -> TruthTrigRates.py truth skim (pathena, SkimCode AthAnalysis 25.2.89) -> rucio download.

Resumable: state in grid/chain_state.json. Polls bigpanda every POLL s. Prints one line per state
change. EXITS (code 0) when every config is downloaded, or (code 1) as soon as any task ends in a
non-clean state (finished-with-failures, failed, broken, aborted, exhausted) or a submission fails,
so a human/agent can decide. Run from evgen_nomupt_acceptance/:  python3 grid/chain.py
(grid/chain_state.json must first hold {cfg: {"evgen": jediTaskID}} from grid/submit_evgen.log)
"""
import json, os, re, subprocess, sys, time, urllib.request, ssl

WD = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SKIM_DIR = os.path.normpath(os.path.join(WD, "../../SkimCode"))
OUT_BASE = os.path.expanduser("~/usatlasdata/pythia_private_nomupt/grid_download")
STATE = os.path.join(WD, "grid/chain_state.json")
POLL = 600
VER = "v1"
PROXY = "/tmp/x509up_u%d" % os.getuid()
ALRB = ("export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase; set --; "
        "source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh > /dev/null; ")
CLEAN_RUNNING = {"registered", "defined", "assigning", "ready", "pending", "scouting", "scouted",
                 "running", "prepared", "throttled", "submitting", "toretry", "toincexec", "rerefine"}
BAD = {"finished", "failed", "broken", "aborted", "exhausted", "tobroken", "toabort"}

# config -> evgen jediTaskID: read from grid/chain_state.json (written after submit_evgen.sh)
EVGEN = {}


def log(msg):
    print(time.strftime("%Y-%m-%d %H:%M:%S"), msg, flush=True)


def panda(task):
    ctx = ssl.create_default_context(capath="/etc/grid-security/certificates")
    ctx.load_cert_chain(PROXY, PROXY)
    req = urllib.request.Request("https://bigpanda.cern.ch/task/%d/?json" % task,
                                 headers={"Accept": "application/json"})
    d = json.load(urllib.request.urlopen(req, context=ctx, timeout=120))
    status = d["task"]["status"]
    outs = [x["datasetname"] for x in d.get("datasets", []) if x.get("type") == "tmpl_output"]
    return status, outs


def run(cmd, cwd):
    p = subprocess.run(["bash", "-c", cmd], cwd=cwd, capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def submit(cmd, cwd, tag):
    rc, out = run(cmd, cwd)
    open(os.path.join(WD, "grid", "submit_%s.log" % tag), "w").write(out)
    m = re.search(r"new jediTaskID=(\d+)", out)
    if rc != 0 or not m:
        log("SUBMIT FAILED %s (rc=%d), see grid/submit_%s.log" % (tag, rc, tag)); sys.exit(1)
    return int(m.group(1))


def phys(cfg):
    jo = [f for f in os.listdir(os.path.join(WD, "run", cfg)) if f.startswith("mc.")][0]
    return jo[3:-3]


def main():
    st = json.load(open(STATE)) if os.path.exists(STATE) else {c: {"evgen": t} for c, t in EVGEN.items()}
    while True:
        for cfg, s in st.items():
            if s.get("downloaded"):
                continue
            stage = "skim" if "skim" in s else "deriv" if "deriv" in s else "evgen"
            status, outs = panda(s[stage])
            if s.get("last_" + stage) != status:
                log("%-26s %-5s task %d -> %s" % (cfg, stage, s[stage], status))
                s["last_" + stage] = status
            if status in BAD:
                log("STOP: %s %s task %d ended '%s'" % (cfg, stage, s[stage], status))
                json.dump(st, open(STATE, "w"), indent=1); sys.exit(1)
            if status != "done":
                if status not in CLEAN_RUNNING:
                    log("note: unrecognised status '%s'" % status)
                continue
            p = phys(cfg)
            if stage == "evgen":
                evnt = [o for o in outs if o.endswith("_EVNT.pool.root/")][0]
                s["evgen_log_ds"] = [o for o in outs if o.endswith("_log.generate/")][0]
                d = os.path.join(WD, "grid/deriv_stage", cfg); os.makedirs(d, exist_ok=True)
                s["deriv"] = submit(ALRB + "asetup Athena,25.0.72 > /dev/null; lsetup panda > /dev/null; "
                    "pathena --trf 'Derivation_tf.py --inputEVNTFile %%IN --outputDAODFile output.pool.root "
                    "--formats TRUTH0' --extOutFile DAOD_TRUTH0.output.pool.root --inDS %s "
                    "--outDS user.yuhang.nomuptAcc.%s.TRUTH0.%s --nFilesPerJob 20" % (evnt, p, VER), d, "deriv_" + cfg)
                log("%-26s deriv submitted: task %d (input %s)" % (cfg, s["deriv"], evnt))
            elif stage == "deriv":
                daod = [o for o in outs if "EXT0" in o][0]
                s["skim"] = submit("source %s/setup_25.sh > /dev/null 2>&1; cd %s/run_evgen_high_pT_additional; "
                    "lsetup panda > /dev/null; pathena TruthTrigRates.py --inDS %s "
                    "--outDS user.yuhang.NTUP.nomuptAcc.%s.%s --nFilesPerJob 20 --mergeOutput"
                    % (SKIM_DIR, SKIM_DIR, daod, p, VER), SKIM_DIR, "skim_" + cfg)
                log("%-26s skim submitted: task %d (input %s)" % (cfg, s["skim"], daod))
            else:
                ntup = [o for o in outs if not o.rstrip("/").endswith(".log")]
                out = os.path.join(OUT_BASE, cfg); os.makedirs(out, exist_ok=True)
                rc, txt = run(ALRB + "lsetup rucio > /dev/null; rucio download --dir %s %s %s"
                              % (out, " ".join(ntup), s["evgen_log_ds"]), out)
                open(os.path.join(out, "rucio_download.log"), "w").write(txt)
                if rc != 0:
                    log("STOP: download failed for %s, see %s/rucio_download.log" % (cfg, out)); sys.exit(1)
                s["downloaded"] = out
                log("%-26s DOWNLOADED -> %s" % (cfg, out))
            json.dump(st, open(STATE, "w"), indent=1)
        json.dump(st, open(STATE, "w"), indent=1)
        if all(s.get("downloaded") for s in st.values()):
            log("ALL_DONE"); return
        time.sleep(POLL)


if __name__ == "__main__":
    main()
