import argparse, os, sys
import subprocess
import datetime, time
import psutil
import multiprocessing as mp  
from comp import long_time

def monitor(func):
    worker_process = mp.Process(target=func, args=(10,))
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
        psutil.virtual_memory()
        time.sleep(0.01)
    logs['exec_time'] = time.time() - p.create_time() #.strftime("%Y-%m-%d %H:%M:%S")

    worker_process.join()
    return logs


if __name__ == "__main__":
    log = monitor(long_time)
    print(log)