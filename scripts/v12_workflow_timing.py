"""Append actual public-workflow timings without changing statistical RNG/tasks."""
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import time


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def checkpoint_rows(path):
    if path is None or not Path(path).exists():return 0
    with Path(path).open('rb') as stream:
        return max(0,sum(line.endswith(b'\n') for line in stream)-1)


class WorkflowTimer:
    """Each line is one actual invocation, including errors and partial resumes.

    Completed checkpoint rows are observations, not a count of optimizer calls.
    Lifetime RSS is explicitly cumulative and cannot be interpreted as the
    exclusive peak memory of a phase. Timings of concurrent workflows overlap.
    """
    def __init__(self,path,**metadata):
        self.path=Path(path)
        self.metadata=dict(metadata,host=platform.node(),platform=platform.platform(),pid=os.getpid(),
            invocation=f'{datetime.now(timezone.utc).isoformat()}:{time.monotonic_ns()}',
            telemetry_source_sha256=sha(__file__),
            timing_scope='API call wall time including internal checkpoint I/O; excludes outer result-table export',
            thread_environment={key:os.environ.get(key) for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']})

    def run(self,phase,call,*,checkpoint=None,**metadata):
        record=dict(self.metadata,phase=phase,**metadata,
            started_utc=datetime.now(timezone.utc).isoformat(),
            checkpoint_rows_before=checkpoint_rows(checkpoint))
        started=time.perf_counter();cpu_started=time.process_time()
        try:
            result=call()
            record['status']='returned'
            for key in ['n_planned','n_failed','n_evals','converged','num_maps','n_trees_planned','n_trees_completed','n_trees_incomplete']:
                value=getattr(result,key,None)
                if value is not None and isinstance(value,(int,float,bool,str)):record[key]=value
            return result
        except BaseException as error:
            record.update(status='raised',error=type(error).__name__+': '+str(error))
            raise
        finally:
            record.update(wall_seconds=time.perf_counter()-started,
                process_cpu_seconds=time.process_time()-cpu_started,
                checkpoint_rows_after=checkpoint_rows(checkpoint))
            try:
                import resource
                record['process_lifetime_peak_rss_bytes']=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*(1 if platform.system()=='Darwin' else 1024)
            except ImportError:record['process_lifetime_peak_rss_bytes']=None
            self.path.parent.mkdir(parents=True,exist_ok=True)
            with self.path.open('a') as stream:
                stream.write(json.dumps(record,sort_keys=True,allow_nan=False)+'\n')
                stream.flush()
