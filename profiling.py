import subprocess
import time
import psutil
import multiprocessing as mp  
import os
import stat

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
        try:
            logs['cpu'].append(p.cpu_percent()) # % cpu usage
            logs['mem_usage'].append(p.memory_full_info().uss)
            logs['mem_percent'].append(p.memory_percent(memtype="rss")) # (total memory usage, percent of memory used)
            # logs['disk_usage'].append(psutil.disk_io_counters())
            psutil.virtual_memory()
            time.sleep(0.01)
        except:
            pass
    logs['exec_time'] = time.time() - p.create_time()

    worker_process.join()
    return logs

def monitor_gzip(file, ori_file):

    # data/headerless\ecoli_100Kb_reads_10x.fasta.headerless.gz 292244
    # data/headerless\ecoli_100Kb_reads_120x.fasta.headerless.gz 3508076
    # data/headerless\ecoli_100Kb_reads_20x.fasta.headerless.gz 586108
    # data/headerless\ecoli_100Kb_reads_40x.fasta.headerless.gz 1170066
    # data/headerless\ecoli_100Kb_reads_5x.fasta.headerless.gz 146953
    # data/headerless\ecoli_100Kb_reads_80x.fasta.headerless.gz 2342520

    unsorted_comp_size =  {
    'data/headerless/ecoli_100Kb_reads_10x.fasta.headerless.gz' : 292244,
    'data/headerless/ecoli_100Kb_reads_120x.fasta.headerless.gz' : 3508076,
    'data/headerless/ecoli_100Kb_reads_20x.fasta.headerless.gz' : 586108,
    'data/headerless/ecoli_100Kb_reads_40x.fasta.headerless.gz' : 1170066,
    'data/headerless/ecoli_100Kb_reads_5x.fasta.headerless.gz' : 146953,
    'data/headerless/ecoli_100Kb_reads_80x.fasta.headerless.gz' : 2342520
    }

    logs = {}
    worker_process = subprocess.Popen(["gzip", "-f", "-k", file])
    p = psutil.Process(worker_process.pid)
    # log cpu usage of `worker_process` every 10 ms
    logs = {
        'cpu': [],
        'mem_usage': [],
        'mem_percent': [],
        'disk_usage': [],
    }
    while worker_process.poll() is None:
        try:
            logs['cpu'].append(p.cpu_percent()) # % cpu usage
            logs['mem_usage'].append(p.memory_full_info().uss)
            logs['mem_percent'].append(p.memory_percent(memtype="rss")) # (total memory usage, percent of memory used)
            # logs['disk_usage'].append(psutil.disk_io_counters())
            psutil.virtual_memory()
            time.sleep(0.01)
        except:
            pass
    logs['exec_time'] = time.time() - p.create_time()
    
    logs['compression_ratio'] =  unsorted_comp_size[ori_file] / os.path.getsize(file + ".gz")
    worker_process.wait()
    # os.remove(f'{file}.gz')

    return logs

def long_time(n):
    for i in range(n):
        for j in range(100000):
            i*j

if __name__ == "__main__":
    print('function logs: ', monitor(long_time(500)))
    print('Gzip logs: ', monitor_gzip('data/ecoli_100Kb_reads_80x.fasta', 'data/headerless/ecoli_100Kb_reads_80x.fasta.headerless.gz'))
