import subprocess
import time
import psutil
import multiprocessing as mp  
from comp import long_time
import os

def monitor(func):
    worker_process = mp.Process(target=func)
    worker_process.start()
    p = psutil.Process(worker_process.pid)

    # log cpu usage of `worker_process` every 10 ms
    logs = {
        'cpu': [],
        'mem_usage': [],
        'mem_percent': [],
        'disk_usage': [],
    }
    while worker_process.is_alive():
        logs['cpu'].append(p.cpu_percent()) # % cpu usage
        logs['mem_usage'].append(p.memory_full_info().uss)
        logs['mem_percent'].append(p.memory_percent(memtype="rss")) # (total memory usage, percent of memory used)
        logs['disk_usage'].append(psutil.disk_io_counters())
        psutil.virtual_memory()
        time.sleep(0.01)
    logs['exec_time'] = time.time() - p.create_time()

    worker_process.join()
    return logs

def monitor_gzip(file):
    logs = {}
    worker_process = subprocess.Popen(["gzip", "-f", "-k", file])
    p = psutil.Process(worker_process.pid)
    print(p)
    # log cpu usage of `worker_process` every 10 ms
    logs = {
        'cpu': [],
        'mem_usage': [],
        'mem_percent': [],
        'disk_usage': [],
    }
    print(worker_process)
    while worker_process.poll() is None:
        logs['cpu'].append(p.cpu_percent()) # % cpu usage
        logs['mem_usage'].append(p.memory_full_info().uss)
        logs['mem_percent'].append(p.memory_percent(memtype="rss")) # (total memory usage, percent of memory used)
        logs['disk_usage'].append(psutil.disk_io_counters())
        psutil.virtual_memory()
        time.sleep(0.01)
    logs['exec_time'] = time.time() - p.create_time()
    worker_process.wait()
    print('Gzip logs: ', logs)
    return logs

if __name__ == "__main__":
    log = monitor(long_time(5))
    print(log)
    monitor_gzip('data/ecoli_100Kb_reads_80x.fasta')
    os.remove('data/ecoli_100Kb_reads_80x.fasta.gz')
