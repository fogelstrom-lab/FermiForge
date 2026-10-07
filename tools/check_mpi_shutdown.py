"""Isolated local-runtime experiment; suspends only its own test launcher.

No solver runs are touched. Compare clean exits and delayed-server shutdown
with two PMIx acknowledgement timeouts. This intentionally bypasses Open
MPI's finalize fence ONLY in the probe so the injected delay reaches the ACK
stage. It demonstrates a failure mechanism, not the cause of production delay.
Output is kept in a unique work folder.
"""
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]


def main():
    folder = Path(tempfile.mkdtemp(prefix="mpi-shutdown-", dir=ROOT / "work"))
    compiler = shutil.which("mpifort") or "/opt/homebrew/bin/mpifort"
    launcher = shutil.which("mpiexec") or "/opt/homebrew/bin/mpiexec"
    exe = folder / "probe"
    subprocess.run([compiler, "-J", str(folder), str(ROOT / "src/mpi_shutdown_audit.f90"),
                    str(ROOT / "tests/mpi_shutdown_probe.f90"), "-o", str(exe)],
                   cwd=folder, check=True)
    results = []
    for label, timeout, delay in [("baseline",10,0), ("delayed-default",10,15),
                                  ("delayed-extended",60,15)]:
        prefix = folder / label
        env = dict(os.environ, OMPI_MCA_pml="ob1", OMPI_MCA_btl="self,sm",
                   # Fault injection only: bypass the internal PMIx fence so
                   # the pause exercises the subsequent finalize ACK timeout.
                   # NEVER apply this setting to production solver runs.
                   OMPI_MCA_async_mpi_finalize="1",
                   PMIX_MCA_pmix_finalize_timeout=str(timeout),
                   PMIX_MCA_pmix_client_base_verbose="2")
        command = [launcher,"--host","localhost:2","--map-by","ppr:2:node",
                   "--bind-to","none","-n","2",str(exe),str(prefix)]
        with prefix.with_suffix(".log").open("w") as log:
            proc = subprocess.Popen(command, env=env, stdout=log, stderr=subprocess.STDOUT)
            stopped = False
            try:
                deadline = time.monotonic()+30
                while not all(Path(str(prefix)+f".ready.{r}").exists() for r in range(2)):
                    if proc.poll() is not None or time.monotonic()>deadline:
                        raise RuntimeError(f"probe did not become ready: {prefix}.log")
                    time.sleep(0.05)
                if delay:
                    os.kill(proc.pid,signal.SIGSTOP)
                    stopped = True
                    time.sleep(delay)
                    os.kill(proc.pid,signal.SIGCONT)
                    stopped = False
                code = proc.wait(timeout=90)
            finally:
                if stopped:
                    os.kill(proc.pid,signal.SIGCONT)
                if proc.poll() is None:
                    proc.terminate()
                    proc.wait(timeout=20)
        text = prefix.with_suffix(".log").read_text()
        audits = [Path(str(prefix)+f".shutdown.rank{r:06d}.txt").read_text()
                  for r in range(2)]
        result = dict(case=label,exit_code=code,
                      injected_async_finalize=True, server_pause_seconds=delay,
                      pmix_finalize_timeout=timeout,
                      timeout_messages=text.count("finwait timeout fired"),
                      improper_exit="exiting improperly" in text or "exiting\nimproperly" in text,
                      audits=audits)
        results.append(result)
        (folder / "results.json").write_text(json.dumps(results,indent=2)+"\n")
        print(label, "exit",code,"timeouts",result["timeout_messages"],flush=True)
    print("Results:",folder,flush=True)


if __name__ == "__main__":
    main()
