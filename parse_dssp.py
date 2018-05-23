# NOTE: This requires the DSSPData class shown in this gist: https://gist.github.com/jlhg/5181883
import os
import pandas as pd
from multiprocessing import Process, Queue
from DSSPData import DSSPData

datapath = "/home/dillon/data/dssp/raw/"
STATUS_FLAG = -100

def _parse(filename):
    '''Parse a single .dssp file and return its properties.
    Returns:
        dssp_id                 - str
        seq                     - str
        secondary structures    - str
        tco                     - list
        kappa                   - list
        alpha                   - list
        phi                     - list
        psi                     - list
    '''

    dssp = DSSPData()
    dssp.parseDSSP(filename)

    dssp_id = os.path.splitext(os.path.basename(filename))[0]
    seq = ''.join(dssp.aa)
    structs = []
    for i in range(len(dssp.aa)):
        s = dssp.struct[i][2]
        if s == ' ':
            s = 'U'
        structs.append(s)

    ss = ''.join(structs)
    tco = ','.join(dssp.getTCO())
    kappa = ','.join(dssp.getKAPPA())
    alpha = ','.join(dssp.getALPHA())
    phi = ','.join(dssp.getPHI())
    psi = ','.join(dssp.getPSI())

    return (dssp_id, seq, ss, tco, kappa, alpha, phi, psi)

def parser(worker_queue, done_queue, id):
    '''Parser task that parses DSSP files from the worker_queue and places
    the results in the done_queue.
    '''

    local_count = 0

    while True:
        filename = worker_queue.get()
        # check if done
        if filename is None:
            return

        local_count += 1
        if local_count % 1000 == 0:
            print("Process #%d has handled %d records." % (id, local_count))

        res = _parse(filename)
        done_queue.put(res)


if __name__ == '__main__':

    num_workers = 4
    files = os.listdir(datapath)
    num_records = len(files)
    num_files = 10
    tenth = num_records // num_files

    print("Processing %d individual records into %d files." % (num_records, num_files))

    worker_queue = Queue()
    done_queue = Queue()

    print("Controller's pid is %d" % (os.getpid()))
    print("Spawning %d workers..." % (num_workers))
    workers = []
    for i in range(num_workers):
        p = Process(target=parser, args=(worker_queue, done_queue, i))
        workers.append(p)
        p.start()

    # fill the queue with files
    for filename in files:
        worker_queue.put(datapath+filename)

    # Add sentinels to signal end
    for _ in range(num_workers):
        worker_queue.put(None)

    print("Parsing DSSP files...")
    results = []
    file_index = 1
    for i in range(num_records):
        res = done_queue.get()
        results.append(res)

        assert len(results) != 0

        # split the data into 10 different files
        if len(results) % tenth == 0:
            print("%2d%% done, saving to dssp_%d.csv" % (file_index*10, file_index))
            filename = "dssp_"+str(file_index)+".csv"
            df = pd.DataFrame.from_records(results, columns=["dssp_id", "seq", "ss", "tco", "kappa", "alpha", "phi", "psi"])
            df.to_csv(filename, index=False)
            file_index += 1
            results = []

    print("Finished. Saving remaining %d to dssp_%d.dssp" % (len(results), file_index))
    filename = "dssp_"+str(file_index)+".csv"
    df = pd.DataFrame.from_records(results, columns=["dssp_id", "seq", "ss", "tco", "kappa", "alpha", "phi", "psi"])
    df.to_csv(filename, index=False)

    print("Joining %d workers..." % (num_workers))
    for p in workers:
        p.join()
