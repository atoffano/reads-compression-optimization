import time
import psutil
import multiprocessing as mp  
import os, gzip
import shutil
import pca_sort

def monitor(func, input_file, compare_to):

    unsorted_comp_size =  {
    'data/headerless/ecoli_100Kb_reads_10x.fasta.headerless.gz' : 292244,
    'data/headerless/ecoli_100Kb_reads_120x.fasta.headerless.gz' : 3508076,
    'data/headerless/ecoli_100Kb_reads_20x.fasta.headerless.gz' : 586108,
    'data/headerless/ecoli_100Kb_reads_40x.fasta.headerless.gz' : 1170066,
    'data/headerless/ecoli_100Kb_reads_5x.fasta.headerless.gz' : 146953,
    'data/headerless/ecoli_100Kb_reads_80x.fasta.headerless.gz' : 2342520
    }

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

    p.wait()
    worker_process.join()

    with open(input_file, 'rb') as f_in:
        with gzip.open(f"{input_file.replace('data/', '')}.gz", 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)    
    logs['compression_ratio'] =  unsorted_comp_size[compare_to] / os.path.getsize(f'{input_file}.gz')

    return logs



if __name__ == "__main__":
    print(monitor(
        func=pca_sort.sort_by_pca('data/ecoli_100Kb_reads_80x.fasta', "out_x.fasta", 50000),
        input_file='out_x.fasta',
        compare_to='data/headerless/ecoli_100Kb_reads_80x.fasta.headerless.gz')
        )
    os.remove('out_x.fasta')
    os.remove('out_x.fasta.gz')
